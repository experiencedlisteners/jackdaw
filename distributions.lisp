(cl:in-package #:jackdaw)

(defun hash-table->alist (hashtab &key (key #'identity) (value #'identity))
  (let ((alist '()))
    (maphash #'(lambda (k v)
		 (push (cons (funcall key k)
			     (funcall value v))
		       alist))
             hashtab)
    alist))

(defun alist->hash-table (alist &key (test #'equal))
  (let ((hashtable (make-hash-table :test test)))
    (mapc #'(lambda (x)
	      (setf (gethash (car x) hashtable) (cdr x)))
          alist)
    hashtable))

(defmacro defdistribution (class superclasses parameters
			   (&optional symbol arguments distribution)
			   &body body)
  "Define a new probability distribution.

Define a class named CLASS with superclasses SUPERCLASSES. One of the superclasses must
be a subclass of PROBABILITY-DISTRIBUTION. If no superclass is given, the class will be created as
a subclass of PROBABILITY-DISTRIBUTION.

A list of parameters can be given in the PARAMETERS argument. 
The probability function should be defined in BODY. When given, SYMBOL, ARGUMENTS, and 
DISTRIBUTION are available within BODY and are bound to respectively the symbol of which 
the probability is required, the values of the variables on which the distribution is 
conditioned, and the distribution instance itself."
  (let* ((direct-slots (%lambda-list->direct-slots parameters)))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c)
					 (subtypep c 'probability-distribution))
				       superclasses))
		    (1- (length superclasses))))
	    ()
	    "At least one superclass of ~a must be of type PROBABILITY-DISTRIBUTION."
	    class)
    `(%muffle-redefinition-warnings
       (defclass ,class ,(if (null superclasses) '(probability-distribution) superclasses)
	 ((%parameters :initform ',(mapcar #'%param-name (remove '&key parameters)))	  
	  ,@direct-slots))
       (defmethod probability ((d ,class) observation)
	 (pr:in (let (,@(unless (null symbol) `((,symbol (car observation))))
		      ,@(unless (null arguments) `((arguments (cdr observation))))
		      ,@(unless (null distribution) `((,distribution d)))
		      ,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
			      collect `(,p (,p d))))
		  (declare (ignorable ,@(mapcar #'%param-name (remove '&key parameters))))
		  ,(if (and (listp arguments) (not (null arguments)))
		       `(destructuring-bind ,arguments
			    arguments
			  ,@body)
		       `(let (,@(unless (null arguments) `((,arguments arguments))))
			  ,@body)))))
       (defun ,(intern (format nil "MAKE-~A-DISTRIBUTION" (symbol-name class))) ,parameters
	 (make-instance ',class ,@(%lambda-list->plist parameters))))))

(defmacro defestimator (distribution (&optional data distribution-symbol)
			(&optional symbol arguments)
			binding-spec parameter-setters
			&key dataset-handler sequence-handler observation-handler)
  "Define an estimator of DISTRIBUTION.

Estimators estimate the parameters of a distribution based on a set of observations."
  (let* ((observation-handler-code
	   `(dolist (observation sequence)
	      (let (,@(unless (null symbol) `((,symbol (car observation))))
		    ,@(unless (null arguments) `((,arguments (cdr observation)))))
		,observation-handler)))
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
    `(defmethod estimate ((d ,distribution) data)
       (let (,@binding-spec
	     ,@(unless (null data) `((,data data)))
	     ,@(unless (null distribution-symbol) `((,distribution-symbol d))))
	 ,@(unless (every #'null (list dataset-handler sequence-handler observation-handler))
	    dataset-handler-code)
	 ,@(loop for (parameter setter) in parameter-setters collect
		 `(setf (slot-value d ',parameter)
			,setter)))
       d)))

(defmethod probability :before ((d probability-distribution) observation)
  (dolist (p (slot-value d '%parameters))
    (unless (slot-boundp d p)
      (error "Distribution parameter ~a of ~a not initialized. You
must either estimate the model or provide the parameters manually before
attempting to access probabilities." p (type-of d)))))

(defmethod conditional-probability ((d probability-distribution) observation
				    &optional arguments)
  (declare (ignore arguments))
  (probability d observation))

(defmethod conditional-probability ((d conditional-probability-distribution) observation
				    &optional arguments)
  (probability d (cons observation arguments)))

(defmethod conditional-probabilities ((d probability-distribution) congruent-values
				      &optional arguments)
  "Obtain the probabilities of a list of possible values given arguments."
  (mapcar (lambda (val) (conditional-probability d val arguments))
	  congruent-values))

(defmethod estimate ((d probability-distribution) data)
  (error "Distribution of type ~a has no defined estimator." (type-of d)))

;;;;;;;;;;; PROBABILITY DISTRIBUTIONS ;;;;;;;;

;;;;;;;;;;;;;;;;;;; BERNOUILLI ;;;;;;;;;;;;;;;

(defdistribution bernouilli () (p &key (psymbol t))
    (symbol)
  (if (equal symbol psymbol) p
      (- 1 p)))

(defestimator bernouilli (data) (symbol arguments)
    ((psymbol (car (car (car data))))
     (total-count 0)
     (count 0))
    ((p (/ count total-count))
     (psymbol psymbol))
  :observation-handler
  (progn
    (incf total-count)
    (when (not (null arguments))
      (warn "Arguments defined for Bernouilli distribution but these are ignored."))
    (when (equal symbol psymbol)
      (incf count))))

;;;;;;;;;;;;;;;;;;;; UNIFORM ;;;;;;;;;;;;;;;;;

(defdistribution uniform (conditional-probability-distribution)
    () () 1) ;; probabilities are normalized automatically
(defestimator uniform () () () ())

;;;;;;;; CONDITIONAL PROBABILITY TABLE ;;;;;;;  

(defdistribution cpt (conditional-probability-distribution)
    (&key domain (cpt (make-hash-table))) (symbol args)
  "A conditional probability table."
  (multiple-value-bind (p found?)
      (gethash (cons symbol args) cpt)
    (unless found?
      (warn "Probability of ~a given ~a not found in conditional probability table." symbol args))
    p))

(defestimator cpt (data) (symbol arguments)
    ((counts (make-hash-table :test 'equal))
     (context-counts (make-hash-table :test 'equal))
     (domain))
    ((cpt
      (let ((cpt (make-hash-table :test 'equal)))
	(maphash (lambda (obs count)
		   (setf (gethash obs cpt)
			 (/ count (gethash (cdr obs) context-counts))))
		 counts)
	cpt))
     (domain domain))
  :observation-handler
  (progn
    (unless (member symbol domain :test #'equal)
      (push symbol domain))
    (setf (gethash (cons symbol arguments) counts)
	  (1+ (gethash (cons symbol arguments) counts 0)))
    (setf (gethash arguments context-counts)
	  (1+ (gethash arguments context-counts 0)))))

(defwriter cpt (m) (hash-table->alist (p m)))
(defreader cpt (m p) (setf (slot-value m 'p) (alist->hash-table p)))

(defmethod initialize-instance :after ((d cpt) &key alist-cpt)
  "Parameters can be supplied as an ALIST: a list with items (PARAM . PROB). 
The context of a parameter is (CDR PROB), and corresponds to a list of states 
corresponding to variables that D is conditioned on. If D is not conditioned 
on anything, the context may be set to NIL. This means that each parameter 
must be a list of length 1 (the CDR of which is NIL)."
  (setf (slot-value d 'cpt) (alist->hash-table alist-cpt))
  (setf (slot-value d 'domain)
	(remove-duplicates
	 (loop for param being each hash-key of (cpt d) collect (car param))
	 :test #'equal))

(defmethod probability-table ((d cpt))
  (let* ((table))
    (maphash #'(lambda (k value)
		 (push (append k (list value)) table))
	     (cpt d))
    table))

;;;;;;;;;;;;;;;;;; NGRAM MODEL ;;;;;;;;;;;;;;;

(defdistribution ngram-model (conditional-probability-distribution)
    (&key (cpt (make-cpt-distribution))) (symbol arguments)
  "Observations of an ngram model are of the form (NGRAM . ARGUMENTS), where 
NGRAM is of the form (Xn Xn-1 ... X0). In the body of this definition, NGRAM is 
bound to SYMBOL and ARGUMENTS to ARGUMENTS."
  (probability cpt (append symbol arguments)))

(defestimator ngram-model (data dist) (observation arguments)
	      ((cpt-dataset))
	      ((cpt (estimate (cpt dist) (print cpt-dataset))))
	      :sequence-handler
	      (push nil cpt-dataset)
	      :observation-handler
	      (push (append observation arguments)
		    (car cpt-dataset)))

;;;;;;;;;;;;;;;;;;; IDYOM PPM ;;;;;;;;;;;;;;;;

(defestimator ppm:ppm (data model) () () ()
	      :dataset-handler
	      (ppm:model-dataset model data :construct? t :predict? nil))

;;;;;;;;;;;;;;;;;;;;;; PPMS ;;;;;;;;;;;;;;;;;;

(defclass ppms (conditional-probability-distribution)
  ((alphabet :reader alphabet :initform nil)
   (escape :initarg :escape :reader escape :initform :c)
   (mixtures :initarg :mixtures :reader mixtures :initform t)
   (update-exclusion :initarg :update-exclusion :reader update-exclusion :initform nil)
   (order-bound :initarg :order-bound :reader order-bound :initform nil)
   (ppms :reader ppms :initform (make-hash-table :test #'equal))
   (locations :accessor locations :initform (make-hash-table :test #'equal)))
  (:documentation "A set of PPM models that can be conditioned on other variables.
Which PPM model is used depends on the values of the variables conditioned on."))
	      
(defmethod probability ((d ppms) observation)
  (warn "Calling PROBABILITY directly is inefficient for PPMS. Use PROBABILITIES instead.")
  (let* ((value (car observation))
	 (arguments (cdr observation))
	 (context (cdr value))
	 (model (get-model d arguments))
	 (location (get-location d model context arguments))
	 (alphabet (mapcar #'car congruent-values))
	 (distribution (ppm::get-distribution d location)))
    (find value distribution :key #'car :test #'equal)))

(defestimator ppms (data) (symbol arguments)
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
	(setf (gethash arguments (ppms d)) model))))))

(defmethod spawn-ppm ((d ppms))
  "Create a PPM model with the parameter settings stored in D."
  (make-instance
   'ppm:ppm :escape (escape d) :order-bound (order-bound d)
	    :mixtures (mixtures d) :update-exclusion (update-exclusion d)
	    :normalise nil ;; avoids double work since normalization is done by us
	    :alphabet (alphabet d)))

(defmethod next-sequence ((d probability-distribution) congruent-states)
  "Called after each training sequence. May be used to update
model state.")

(defmethod next-sequence ((d ppms) congruent-values)
  (loop for value in congruent-values do
       (let* ((sequence (car value))
	      (arguments (cdr value))
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
	      
(defmethod conditional-probabilities ((d ppms) possible-values
				      &optional arguments)
  "Obtain the location object of the appropriate PPM model given context.
Context is obtained by accessing the previous self of the current variable, 
which if the current variable is an accumulator, must represent the context.
Note that PARENTS-STATE represents a state in the current moment in which any
parent variables are instantiated."
  (let* ((context (cdr (car possible-values)))
	 (model (get-model d arguments))
	 (location (get-location d model context arguments))
	 (alphabet (mapcar #'car possible-values)))
    (ppm:set-alphabet model alphabet)
    (mapcar (lambda (item) (pr:in (cadr item)))
	    (ppm::get-distribution model location))))
