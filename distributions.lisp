(cl:in-package #:jackdaw)

(defun hash-table->alist (hashtab &key (key #'identity))
  (let ((alist '()))
    (maphash #'(lambda (k value)
		 (setq alist
		       (cons (list k (funcall key value)) alist)))
             hashtab)
    alist))

(defun alist->hash-table (alist &key (test #'equal))
  (let ((hashtable (make-hash-table :test test)))
    (mapc #'(lambda (x) (setf (gethash (car x) hashtable) (cadr x)))
          alist)
    hashtable))

(defmacro defdistribution (class superclasses parameters (args symbol) &body body)
  (let* ((direct-slots (%lambda-list->direct-slots parameters)))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c) (subtypep c 'distribution)) superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type DISTRIBUTION." class)
    `(%muffle-redefinition-warnings
       (defclass ,class ,(if (null superclasses) '(distribution) superclasses)
	 ((%parameter-slots :initform ',(mapcar #'%param-name (remove '&key parameters)))	  
	  ,@direct-slots))
       (defmethod probability ((d ,class) observation)
	 (pr:in (let ((,symbol (car observation))
		      (arguments (cdr observation))
		      ,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
			      collect `(,p (,p d))))
		  (declare (ignorable ,symbol))
		  ,(if (listp args)
		       `(destructuring-bind ,args
			    arguments
			  ,@body)
		       `(let ((,args arguments))
			  (declare (ignorable ,args))
			 ,@body)))))
       (defun ,(intern (format nil "MAKE-~A-DISTRIBUTION" (symbol-name class))) ,parameters
	 (make-instance ',class ,@(%lambda-list->plist parameters))))))

(defmacro defestimator (distribution data symbol arguments
			binding-spec parameter-setters
			&key dataset-handler sequence-handler observation-handler)
  (let* ((observation-handler-code
	   `(dolist (observation sequence)
	      (let ((,symbol (car observation))
		    (,arguments (cdr observation)))
		,@observation-handler)))
	 (sequence-handler-code
	   `(dolist (sequence ,data)
	      ,@(unless (null sequence-handler) (list sequence-handler))
	      ,@(unless (null observation-handler)
		 (list observation-handler-code))))
	 (dataset-handler-code
	   `(,@(unless (null dataset-handler) (list dataset-handler))
	     ,@(unless (and (null sequence-handler)
			   (null observation-handler))
		(list sequence-handler-code)))))
    `(defmethod estimate ((d ,distribution) ,data)
       (let ,binding-spec
	 ,@(unless (every #'null (list dataset-handler sequence-handler observation-handler))
	    dataset-handler-code)
	 ,@(loop for (parameter setter) in parameter-setters collect
		 `(setf (slot-value d ',parameter)
			,setter)))
       d)))

(defclass distribution ()
  ((arguments :initarg :arguments :reader arguments :initform nil)
   (variable-symbol :initarg :variable-symbol :reader variable-symbol)
   (%parameters :initform nil)))

(defmethod probability :before ((d distribution) observation)
  (dolist (p (slot-value d '%parameters))
    (unless (slot-boundp d p)
      (error "Distribution parameter ~a of ~a not initialized. You
must either estimate the model or provide the parameters manually before
attempting to access probabilities." p (type-of d)))))

(defmethod estimate ((d distribution) data)
  (error "Distribution of type ~a has no defined estimator." (type-of d)))

;;;;;;;;;;;;;;;;;;; Probability distributions ;;;;;;;;;;;;;;;;;;;

(defdistribution bernouilli () (p &key (psymbol t))
    (args symbol)
  (if (equal symbol psymbol) p
      (- 1 p)))

(defestimator bernouilli data symbol arguments
    ((psymbol (car (car (car data))))
     (total-count 0)
     (count 0))
    ((p (/ count total-count))
     (psymbol psymbol))
  :observation-handler
  ((incf total-count)
   (when (not (null arguments))
     (warn "Arguments defined for Bernouilli distribution but these are ignored."))
   (when (equal symbol psymbol)
     (incf count))))

(defdistribution uniform () () (args symbol) 1) ;; probabilities are normalized automatically
(defestimator uniform data s a nil nil)

(defdistribution deterministic () () (() symbol) 1)
(defestimator deterministic data s a nil nil)
(defmethod probabilities ((d deterministic) parents-state  congruent-states)
  (assert (eq (length congruent-states) 1) ()
	  "Variable with deterministic distribution must have exactly one congruent state.")
  (call-next-method))

(defdistribution cpt () (&key (cpt (make-hash-table))) (args symbol)
  (multiple-value-bind (p found?)
      (gethash (cons symbol args) cpt)
    ;;(format t "P~w = ~a~%" (cons symbol args) (gethash (cons symbol args) cpt))
    (unless found?
      (warn "Probability of ~a given ~a not found in conditional probability table." symbol args))
    p))

(defestimator cpt data symbol arguments
    ((counts (make-hash-table :test 'equal))
     (context-counts (make-hash-table :test 'equal)))
    ((cpt
      (let ((cpt (make-hash-table :test 'equal)))
	(maphash (lambda (obs count)
		   (setf (gethash obs cpt)
			 (/ count (gethash (cdr obs) context-counts))))
		 counts)
	cpt)))
  :observation-handler
  ((setf (gethash (cons symbol arguments) counts)
	 (1+ (gethash (cons symbol arguments) counts 0)))
   (setf (gethash arguments context-counts)
	 (1+ (gethash arguments context-counts 0)))))

(defmethod estimate ((model ppm::ppm) data)
  (ppm:model-dataset model data :construct? t :predict? nil)
  model)

(defclass ppms (distribution)
  ((alphabet :reader alphabet :initform nil)
   (escape :initarg :escape :reader escape :initform :c)
   (mixtures :initarg :mixtures :reader mixtures :initform t)
   (update-exclusion :initarg :update-exclusion :reader update-exclusion :initform nil)
   (order-bound :initarg :order-bound :reader order-bound :initform nil)
   (ppms :reader ppms :initform (make-hash-table :test #'equal))
   (locations :accessor locations :initform (make-hash-table :test #'equal)))
  (:documentation "A set of PPM models that can be conditioned on other variables.
Which PPM model is used depends on the values of the variables conditioned on."))

(defestimator ppms data symbol arguments
  ((ppms (make-hash-table :test #'equal))
   (datasets (make-hash-table :test #'equal)))
  ((ppms (progn
	   (maphash (lambda (context data)
		      (let ((ppm (spawn-ppm d)))
			(estimate ppm data)
			(setf (gethash context ppms) ppm)))
		    datasets)
	  ppms)))
  :sequence-handler
  (let ((context (cdr (car sequence)))
	(observation-sequence))
    (dolist (observation sequence)
      (when (not (equal (cdr observation) context))
	(error "Context change within a sequences not allowed for a sequence model."))
      (unless (inactive? (car observation)) (push (caar observation) observation-sequence)))
    ;;(print (coerce (reverse observation-sequence) 'string))
    (when (null	(gethash context datasets))
      (setf (gethash context datasets) nil))
    (push (reverse observation-sequence) (gethash context datasets))))

(defwriter distribution (m)
	   (loop for s in (%parameter-slots m) collect (slot-value m s)))		 
(defreader distribution (m data)
  (loop for v in data for s in (%parameter-slots m)
	collect (slot-value m s)))
(defwriter cpt (m) (hash-table->alist (p m)))
(defreader cpt (m p) (setf (slot-value m 'p) (alist->hash-table p)))
(defwriter ppm:ppm (m)
    (list :leaves (utils:hash-table->alist (ppm::ppm-leaves m))
	  :branches (utils:hash-table->alist (ppm::ppm-branches m))
	  :dataset (ppm::dataset->alist m)
	  :alphabet (ppm::ppm-alphabet m)
	  :order-bound (ppm::ppm-order-bound m)
	  :mixtures (ppm::ppm-mixtures m)
	  :escape (ppm::ppm-escape m)
	  :update-exclusion (ppm::ppm-update-exclusion m)))
(defreader ppm:ppm (m data)
  (setf (slot-value m 'ppm::leaves)
	(utils:alist->hash-table (getf data :leaves))
	(slot-value m 'ppm::branches)
	(utils:alist->hash-table (getf data :branches))
	(slot-value m 'ppm::dataset)
	(ppm::alist->dataset (getf data :dataset))
	(slot-value m 'ppm::alphabet)
	(getf data :alphabet)
	(slot-value m 'ppm::order-bound)
	(getf data :order-bound)
	(slot-value m 'ppm::mixtures)
	(getf data :mixtures)
	(slot-value m 'ppm::escape)
	(getf data :escape)
	(slot-value m 'ppm::update-exclusion)
	(getf data :update-exclusion))
    (ppm::initialize m))
(defwriter ppms (d)
    (loop for arguments being the hash-keys of (ppms d) collect
	 (cons arguments (serialize (gethash arguments (ppms d))))))
(defreader ppms (d data)
  (dolist (pair data) ;; Pair structure: (cons arguments serialized-ppm)
    (with-input-from-string (s (write-to-string (cdr pair)))
      (let ((model (make-instance 'ppm:ppm))
	    (arguments (car pair)))
	(deserialize model s)
	;; Store a ppm
	(setf (gethash arguments (ppms d)) model)))))

(defmethod initialize-instance :after ((d cpt) &key alist-cpt)
  "Parameters can be supplied as an ALIST: a list with items (PARAM . PROB). 
The context of a parameter is (CDR PROB), and corresponds to a list of states 
corresponding to variables that D is conditioned on. If D is not conditioned 
on anything, the context may be set to NIL. This means that each parameter 
must be a list of length 1 (the CDR of which is NIL)."
  (dolist (parameter alist-cpt)
    (setf (gethash (car parameter) (cpt d)) (cdr parameter)))
  (let ((contexts (remove-duplicates (mapcar #'cdar alist-cpt) :test #'equal)))
    (dolist (context contexts)
      (let ((sum))
	(maphash (lambda (param v)
		   (when (equal (cdr param) context)
		     (setf sum (apply #'+ (cons v (unless (null sum) (list sum)))))))
		 (cpt d))
	(when (> (abs (- sum 1)) 1.0e-10) ;; Check that sum is approximately one.
	  (warn "Parameters of ~A sum to ~A, not to approximately 1.0, for context ~A."
		(variable-symbol d) sum context))))))

(defmethod initialize-instance :after ((d ppms) &key)
  (loop for p in (arguments d) if (horizontal? p) do
       (warn "~a has previous-moment arguments which is not supported at the
moment. See NEXT-SEQUENCE for PPMS and TRANSITION, which ROTATEs
states." (type-of d))))

(defmethod spawn-ppm ((d ppms))
  (make-instance
   'ppm:ppm :escape (escape d) :order-bound (order-bound d)
	    :mixtures (mixtures d) :update-exclusion (update-exclusion d)
	    :normalise nil ;; since normalization is done by us
	    :alphabet (alphabet d)))

(defmethod next-sequence ((d distribution) congruent-states)
  "Called after each training sequence. May be used to update
model state.")

(defmethod next-sequence ((d ppms) congruent-states)
  (loop for state in congruent-states do
       (let* ((sequence (gethash (variable-symbol d) state))
	      (arguments (mapcar (lambda (v) (gethash v state)) (arguments d)))
	      (model (get-model d arguments)))
	 (update-location d model (cdr sequence) arguments (car sequence))
	 ;; TODO: Is it safe to model the sentinel multiple times for a single model?
	 (update-location d model sequence arguments ppm::*sentinel*)))
  (setf (locations d) (make-hash-table :test #'equal)) ; wipe out locations table
  (loop for k being the hash-keys of (ppms d) do
       (let ((model (gethash k (ppms d))))
	 ;;(when (training? d) (ppm:initialise-virtual-nodes model))
	 (ppm:increment-sequence-front model))))

(defmethod update-location ((d ppms) (model ppm:ppm) context arguments symbol)
  "Find the location associated with CONTEXT and ARGUMENTS (stored 
at (CONS CONTEXT ARGUMENTS)) and update it with SYMBOL. 
Store result at (CONS (CONS SYMBOL CONTEXT) ARGUMENTS)"
  ;;(format t "Updating location at ~a with ~s.~%" context symbol)
  (ppm::add-event-to-model-dataset model symbol)
  (let* ((previous-location (gethash (cons context arguments) (locations d)
				     (ppm:get-root)))
	 (location (ppm::ukkstep model nil previous-location symbol nil))) ;;(training? d))))
    ;;(when (training? d)
      ;;(warn "Updating model with ~a given context ~a" symbol context)
    ;;  (ppm::increment-counts model location))
    (unless (eq symbol ppm::*sentinel*)
      (ppm:increment-event-front model))
    (setf (gethash (cons (cons symbol context) arguments) (locations d)) location)))


(defmethod get-location ((d ppms) (model ppm:ppm) context arguments)
  "Obtain PPM location corresponding to the current context. If not found,
this means that either "
  (multiple-value-bind (location found?)
      (gethash (cons context arguments) (locations d))
    (if found? location
	(if (null context)
	    (ppm:get-root)
	    (update-location d model (cdr context) arguments (car context))))))

(defmethod get-model ((d ppms) arguments)
  "Obtain a PPM model for the current arguments."
  (multiple-value-bind (model found?)
      (gethash arguments (ppms d))
    (unless found?
      (error "No PPM model found for arguments ~a." arguments))
    model))
	      
(defmethod probability-distribution ((d ppms) arguments congruent-values)
  "Obtain the location object of the appropriate PPM model given context.
Context is obtained by accessing the previous self of the current variable, 
which if the current variable is an accumulator, must represent the context.
Note that PARENTS-STATE represents a state in the current moment in which any
parent variables are instantiated."
  (let* ((context (cdr (car congruent-values)))
	 (model (get-model d arguments))
	 (location (get-location d model context arguments))
	 (alphabet (mapcar #'car congruent-values)))
    (ppm:set-alphabet model alphabet)
    (mapcar (lambda (item) (pr:in (cadr item)))
	    (ppm::get-distribution model location))))

(defmethod probability-table ((d cpt))
  (let* ((header (append (arguments d) (list (variable-symbol d) 'probability)))
	 (table))
    (maphash #'(lambda (k value)
		 (push (append k (list value)) table))
	     (cpt d))
    (cons header table)))
