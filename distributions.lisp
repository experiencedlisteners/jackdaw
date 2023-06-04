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
  "Define a new probability distribution, a PROBABILITY method specialized on this distribution,
and a factory, MAKE-<CLASS>-distribution to initialize the distribution.

Define a class named CLASS with superclasses SUPERCLASSES. One of the superclasses must
be a subclass of PROBABILITY-DISTRIBUTION. If no superclass is given, the class will be created as
a subclass of PROBABILITY-DISTRIBUTION.

PARAMETERS is a list of distribution parameters. They can be expressed in lambda-list-like
style and will be translated into class slots. See %LAMBDA-LIST->DIRECT-SLOTS for 
more information.

Furthermore, in the PROBABILITY method (see below), the slot-values of these parameters
are bound to the symbols listed in PARAMETERS.

SYMBOL, ARGUMENTS, and DISTRIBUTION are symbols which, when provided are bound in BODY
respectively to the symbol of which we want to know the probability, the arguments on
which it is conditioned, and a the distribution instance itself.

The PROBABILITY method should, given a symbol, arguments, and probability distribution 
return the symbol's conditional probability. It's contents are determined by the BODY
argument.

Example:

The following defines a Bernouilli distribution with two slots: P and PSYMBOL.
P represents the one parameter of the distribution. PSYMBOL is an optional keyword
parameter that defaults to T and indicates which symbol to assign the probability P to.

(defdistribution bernouilli 
  ()                           ; superclasses
  (p &key (psymbol t))         ; parameters
  (obs)                        ; symbol to bind to in probability method
  (if (equal obs psymbol) p    ; body of the probability method
      (- 1 p)))

qAs you can see, the body of the probability method can access OBS, as well as the
parameters P and PSYMBOL.

The example usage below illustrates the usage of MAKE-BERNOUILLI-DISTRIBUTION, and how PSYMBOL
is a keyword parameter of this factory.

> (probability (make-bernouilli-distribution 0.7) t)
0.7
> (probability (make-bernouilli-distribution 0.7 :psymbol 'a) t)
0.3
> (probability (make-bernouilli-distribution 0.7 :psymbol 'a) t)
0.3
"
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
  "Define an ESTIMATE for method DISTRIBUTION.

DISTRIBUTION is a symbol referring to a subclass of PROBABILITY-DISTRIBUTION.

BINDING-SPEC is a list of symbol bindings that are bound before estimation.

PARAMETER-SETTERS is a list of pairs (lists with two items), the left-hand item represents
a slot of DISTRIBUTION, and the right-hand side is a form to which the slot value will
be set.

Estimators can by providing either one, or a combination of DATASET-HANDLER, SEQUENCE-HANDLER,
or OBSERVATION-HANDLER. Each handler will be called respectively on the entire dataset, each sequence or each observation in each sequence.

A dataset is a set of sequences of observations. See [datasets](concepts#datasets) for more information.

The order of operations is as follows: first the symbols in BINDING-SPEC are bound, next the different handlers are called, and finally the PARAMETER-SETTERS are executed to assign the resulting parameters to slots of the distribution.

Example:

Below is an example of an estimator for this the Bernouilli distribution in the example of
DEFDISTRIBUTION.

(defestimator bernouilli 
  (data)                              ; Symbols to bind in ESTIMATE method
  (symbol arguments)                  ; Symbols to bind in observation handler
    ((psymbol (car (car (car data)))) ; Binding-spec for ESTIMATE method
     (total-count 0)                  ; used in this case to initialize some
     (count 0))                       ; parameters used in estimation
  ((p (/ count total-count))          ; Parameter slot setters to be called
   (psymbol psymbol))                 ; after they have been estimated
  :observation-handler                ; Code for updating parameters used for estimation
  (progn
    (incf total-count)
    (when (not (null arguments))
      (warn \"Arguments defined for Bernouilli distribution but these are ignored.\"))
    (when (equal symbol psymbol)
      (incf count))))

The code below illustrates how estimate can be used with a small dataset - consisting of
three sequences with a variable number of osbervations - to estimate a Bernouilli distribution.
Note that we cannot use MAKE-BERNOUILLI-DISTRIBUTION here since that factory requires us to
provide the parameter P.

> (let ((dataset '(((nil t t t) (t nil) (nil t nil))))
        (d (jd:estimate (make-instance 'jd:bernouilli) dataset)))
    (jd:probability d nil))
2/3
"
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

(defmethod estimate ((d probability-distribution) data)
  (error "Distribution of type ~a has no defined estimator." (type-of d)))

;;;;;;;;;;;;;;;;;;; Probability distributions ;;;;;;;;;;;;;;;;;;;

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

(defdistribution uniform () () () 1) ;; probabilities are normalized automatically
(defestimator uniform () () () ())

(defdistribution ngram-model () (&key (cpt (make-cpt-distribution))) (symbol)
  (probability cpt symbol))

(defestimator ngram-model (data dist) (observation)
	      ((cpt-dataset))
	      ((cpt (estimate (cpt dist) (print cpt-dataset))))
	      :sequence-handler
	      (push nil cpt-dataset)
	      :observation-handler
	      (push observation (car cpt-dataset)))

(defdistribution cpt () (&key domain (cpt (make-hash-table))) (symbol args)
  (multiple-value-bind (p found?)
      (gethash (cons symbol args) cpt)
    ;;(format t "P~w = ~a~%" (cons symbol args) (gethash (cons symbol args) cpt))
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

(defestimator ppm:ppm (data model) () () ()
	      :dataset-handler
	      (ppm:model-dataset model data :construct? t :predict? nil))

(defclass ppms (probability-distribution)
  ((alphabet :reader alphabet :initform nil)
   (escape :initarg :escape :reader escape :initform :c)
   (mixtures :initarg :mixtures :reader mixtures :initform t)
   (update-exclusion :initarg :update-exclusion :reader update-exclusion :initform nil)
   (order-bound :initarg :order-bound :reader order-bound :initform nil)
   (ppms :reader ppms :initform (make-hash-table :test #'equal))
   (locations :accessor locations :initform (make-hash-table :test #'equal)))
  (:documentation "A set of PPM models that can be conditioned on other variables.
Which PPM model is used depends on the values of the variables conditioned on."))

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

(defwriter probability-distribution (m)
  (loop for s in (%parameters m) collect (slot-value m s)))
(defreader probability-distribution (m data)
  (loop for v in data for s in (%parameters m)
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
  (setf (slot-value d 'cpt) (alist->hash-table alist-cpt))
  (setf (slot-value d 'domain)
	(remove-duplicates
	 (loop for param being each hash-key of (cpt d) collect (car param))
	 :test #'equal)))

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

(defmethod probabilities ((d probability-distribution) arguments congruent-values)
  "Obtain the probabilities of a list of congruent values given arguments."
  (mapcar (lambda (s) (probability d (cons s arguments)))
	  congruent-values))
	      
(defmethod probabilities ((d ppms) arguments congruent-values)
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
  (let* ((table))
    (maphash #'(lambda (k value)
		 (push (append k (list value)) table))
	     (cpt d))
    table))
