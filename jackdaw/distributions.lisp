(cl:in-package #:jackdaw)

(defmacro defdistribution (class superclasses parameters (args symbol) &body body)
  (let* ((direct-slots (%lambda-list->direct-slots parameters class)))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c) (subtypep c 'distribution)) superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type DISTRIBUTION." class)
    `(muffle-redefinition-warnings
       (defclass ,class ,(if (null superclasses) '(distribution) superclasses)
	 ((%parameter-slots :initform ',(mapcar #'%param-name (remove '&key parameters)))	  
	  ,@direct-slots))
       (defmethod probability ((d ,class) observation)
	 (pr:in (let ((,symbol (car observation))
		      (arguments (cdr observation))
		      ,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
			      collect `(,p (,p d))))
		  ,(if (listp args)
		       `(destructuring-bind ,args
			    arguments
			  ,@body)
		       `(let ((,args arguments))
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

(defmethod estimate ((d distribution) data)
  (error "Distribution of type ~a has no defined estimator." (type-of d)))

;;;;;;;;;;;;;;;;;;; Probability distributions ;;;;;;;;;;;;;;;;;;;

(defdistribution bernouilli () (p &key (psymbol t))
    (args symbol)
  (let ((p (if (null args) p
	       (cdr (assoc args p :test #'equal)))))
    (if (equal symbol psymbol) p
	(- 1 p))))

(defestimator bernouilli data symbol arguments
    ((symbol (car (car data)))
     (total-count 0)
     (count 0))
    ((p (/ count total-count))
     (psymbol symbol))
  :observation-handler
  ((incf total-count)
   (when (not (null arguments))
     (warn "Arguments defined for Bernouilli distribution but these are ignored."))
   (when (equal symbol symbol)
     (incf count))))

(defdistribution uniform () () (() symbol) 1) ;; probabilities are normalized automatically
(defestimator uniform data s a nil nil)

(defdistribution deterministic () () (() symbol) 1)
(defestimator deterministic data s a nil nil)
(defmethod probabilities ((d deterministic) parents-state  congruent-states)
  (assert (eq (length congruent-states) 1) ()
	  "Variable with deterministic distribution must have exactly one congruent state.")
  (call-next-method))

(defdistribution cpt () (cpt) (args symbol)
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

(defclass accumulator-model (distribution)
  ((alphabet :reader alphabet :initform nil)
   (escape :reader escape :initform :c)
   (mixtures :reader mixtures :initform t)
   (update-exclusion :reader update-exclusion :initform nil)
   (order-bound :reader order-bound :initform nil)
   (ppms :reader ppms :initform (make-hash-table :test #'equal))
   (locations :accessor locations :initform (make-hash-table :test #'equal))
   (training? :accessor training? :initarg :training? :initform nil))
  (:documentation "PPM model that's also a jackdaw-native
 sequence model. The overridden fields serve to set some defaults."))

(defestimator accumulator-model data symbol arguments
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
      (push (caar observation) observation-sequence))
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
(defwriter accumulator-model (d)
    (loop for arguments being the hash-keys of (ppms d) collect
	 (cons arguments (serialize (gethash arguments (ppms d))))))
(defreader accumulator-model (d data)
  (dolist (pair data) ;; Pair structure: (cons arguments serialized-ppm)
    (with-input-from-string (s (write-to-string (cdr pair)))
      (let ((model (make-instance 'ppm:ppm))
	    (arguments (car pair)))
	(deserialize model s)
	;; Store a ppm
	(setf (gethash arguments (ppms d)) model)))))

(defmethod initialize-instance :after ((d cpt) &key parameters)
  "Parameters can be supplied as an ALIST: a list with items (PARAM . PROB). 
The context of a parameter is (CDR PROB), and corresponds to a list of states 
corresponding to variables that D is conditioned on. If D is not conditioned 
on anything, the context may be set to NIL. This means that each parameter 
must be a list of length 1 (the CDR of which is NIL)."
  (dolist (parameter parameters)
    (setf (gethash (car parameter) (p d)) (cdr parameter)))
  (let ((contexts (remove-duplicates (mapcar #'cdar parameters) :test #'equal)))
    (dolist (context contexts)
      (let ((sum))
	(maphash (lambda (param v)
		   (when (equal (cdr param) context)
		     (setf sum (apply #'+ (cons v (unless (null sum) (list sum)))))))
		 (p d))
	(when (> (abs (- sum 1)) 1.0e-10) ;; Check that sum is approximately one.
	  (warn "Parameters of ~A sum to ~A, not to approximately 1.0, for context ~A."
		(dist-var d) sum context))))))

(defmethod initialize-instance :after ((d accumulator-model) &key)
  (loop for p in (arguments d) if (previous? p) do
       (warn "~a has previous-moment arguments which is not supported at the
moment. See NEXT-SEQUENCE for ACCUMULATOR-MODEL and TRANSITION, which ROTATEs
states." (type-of d))))

(defmethod spawn-ppm ((d accumulator-model))
  (make-instance
   'ppm:ppm :escape (escape d) :order-bound (order-bound d)
   :mixtures (mixtures d) :update-exclusion (update-exclusion d)
   :alphabet (alphabet d)))

(defmethod next-sequence ((d distribution) congruent-states)
  "Called after each training sequence. May be used to update
model state.")

(defmethod next-sequence ((d accumulator-model) congruent-states)
  (loop for state in congruent-states do
       (let* ((sequence (getarg (previous (dist-var d)) state))
	      (arguments (mapcar (lambda (v) (getarg (previous v) state))
				 (remove (dist-var d) (arguments d))))
	      (model (get-model d arguments)))
	 (update-location d model (cdr sequence) arguments (car sequence))
	 ;; TODO: Is it safe to model the sentinel multiple times for a single model?
	 (update-location d model sequence arguments ppm::*sentinel*)))
  (setf (locations d) (make-hash-table :test #'equal)) ; wipe out locations table
  (loop for k being the hash-keys of (ppms d) do
       (let ((model (gethash k (ppms d))))
	 (when (training? d) (ppm:initialise-virtual-nodes model))
	 (ppm:increment-sequence-front model))))

(defmethod update-location ((d accumulator-model) (model ppm:ppm) context arguments symbol)
  "Find the location associated with CONTEXT and ARGUMENTS (stored 
at (CONS CONTEXT ARGUMENTS)) and update it with SYMBOL. 
Store result at (CONS (CONS SYMBOL CONTEXT) ARGUMENTS)"
  ;;(format t "Updating location at ~a with ~s.~%" context symbol)
  (ppm::add-event-to-model-dataset model symbol)
  (let* ((previous-location (gethash (cons context arguments) (locations d)
				     (ppm:get-root)))
	 (location (ppm::ukkstep model nil previous-location symbol (training? d))))
    (when (training? d)
      ;;(warn "Updating model with ~a given context ~a" symbol context)
      (ppm::increment-counts model location))
    (unless (eq symbol ppm::*sentinel*)
      (ppm:increment-event-front model))
    (setf (gethash (cons (cons symbol context) arguments) (locations d)) location)))


(defmethod get-location ((d accumulator-model) (model ppm:ppm) context arguments)
  "Obtain PPM location corresponding to the current context. If not found,
this means that either "
  (multiple-value-bind (location found?)
      (gethash (cons context arguments) (locations d))
    (if found? location
	(if (null context)
	    (ppm:get-root)
	    (update-location d model (cdr context) arguments (car context))))))

(defmethod get-model ((d accumulator-model) arguments)
  "Obtain a PPM model for the current arguments. If it doesn't exist,
create one and a corresponding root location for the empty context ()."
  (multiple-value-bind (model found?)
      (gethash arguments (ppms d))
    (if found? model
	(let ((model (spawn-ppm d)))
	  (setf (gethash arguments (ppms d)) (spawn-ppm d))
	  model))))

(defmethod get-distribution ((d accumulator-model) table arguments congruent-states)
  "Obtain the location object of the appropriate PPM model given context.
Context is obtained by accessing the previous self of the current variable, 
which if the current variable is an accumulator, must represent the context.
Note that PARENTS-STATE represents a state in the current moment in which any
parent variables are instantiated."
  (let* ((context (cdr (car congruent-states)))
	 (model (get-model d arguments))
	 (location (get-location d model context arguments))
	 (alphabet (mapcar #'car congruent-states)))
    (ppm:set-alphabet model alphabet)
    (dolist (p (ppm::get-distribution model location))
      (let ((symbol (car p))
	    (probability (pr:in (cadr p))))
	(cons symbol context)
	(setf (gethash (cons symbol context) table) probability)))
    table))

(defmethod probabilities ((d distribution) parents-state congruent-states)
  "Obtain the probability of provided CONGRUENT-STATES by the PROBABILITY method.
This is just a wrapper for GET-DISTRIBUTION which grabs arguments from PARENTS-STATE
and avoids a call to GET-DISTRIBUTION when the variable is inactive."
  (let ((table (make-hash-table :test #'equal))
	(arguments (mapcar (lambda (v) (getarg v parents-state)) (arguments d))))
    (if (equal congruent-states (list +inactive+))
	(setf (gethash +inactive+ table) (pr:in 1))
	(get-distribution d table arguments congruent-states))
    table))

(defmethod get-distribution ((d distribution) table arguments congruent-states)
  (let* ((probabilities (mapcar (lambda (s) (probability d (cons s arguments)))
				congruent-states))
	 (sum (apply #'pr:add probabilities)))
    (loop for s in congruent-states for p in probabilities do
	 (setf (gethash s table) (pr:div p sum)))
    table))
