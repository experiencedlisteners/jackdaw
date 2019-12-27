(cl:in-package #:jackdaw)

;;;;;;;;;;;;;;;;;;; Probability distributions ;;;;;;;;;;;;;;;;;;;

(defclass distribution ()
  ((arguments :initarg :arguments :reader arguments :initform nil)
   (variable :initarg :variable :reader dist-var)))
(defclass bernouilli (distribution)
  ((p :initarg :p :accessor p :initform .5)
   (symbols :initarg :symbols :accessor symbols)))
(defclass categorical (distribution)
  ((category-counts :accessor category-counts :initform (make-hash-table :test #'equal))
   (p :reader p :initform (make-hash-table :test #'equal))))
(defclass uniform (distribution) ())
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


(defwriter distribution (m) nil)
(defreader distribution (m d) (declare (ignorable d)))
(defwriter bernouilli (m) (list (p m) (symbols m)))
(defreader bernouilli (m d)
  (setf (slot-value m 'p) (first d)
	(slot-value m 'symbols) (second d)))
(defwriter categorical (m) (hash-table->alist (p m)))
(defreader categorical (m p) (setf (slot-value m 'p) (alist->hash-table p)))
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

(defmethod initialize-instance :after ((d categorical) &key parameters)
  "Parameters must be supplied as an ALIST: a list with items (PARAM . PROB). 
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
  (loop for p in (arguments d) if (apriori? p) do
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
       (let* ((sequence (getarg (apriori (dist-var d)) state))
	      (arguments (mapcar (lambda (v) (getarg (apriori v) state))
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

(defmethod observe ((d distribution) state symbol training?)
  "Ignore observations by default.")

(defmethod observe ((d categorical) state symbol training?)
  (when training?
    (let* ((arguments (mapcar (lambda (v) (getarg v state)) (arguments d)))
	   (counts (category-counts d))
	   (new-arg-count (1+ (gethash arguments counts 0)))
	   (new-s-count (1+ (gethash (cons symbol arguments) counts 0))))
      (setf (gethash arguments counts) new-arg-count)
      (setf (gethash (cons symbol arguments) counts) new-s-count)
      (setf (gethash (cons symbol arguments) (p d))
	    (pr:in (/ new-s-count new-arg-count))))))

(defmethod probability ((d uniform) arguments symbol) (pr:in 1)) ; is normalised later

(defmethod probability ((d bernouilli) arguments symbol)
  (when (not (null arguments))
    (warn "It looks like you're conditioning a Bernouilli distribution on something,
the implementation does not support this."))
  (pr:in (if (equal symbol (car (symbols d)))
	     (p d)
	     (progn
	       (unless (equal symbol (cadr (symbols d)))
		 (warn "Generating (1 - p) probability for unfamiliar symbol."))
	       (- 1 (p d))))))

(defmethod set-param ((d categorical) arguments symbol probability)
  "Set probability of distribution parameter. Probabilities must be given not in
log representation. Caller must ensure probabilities sum to one."
  (setf (gethash (cons symbol arguments) (p d)) probability))

(defmethod probability ((d categorical) arguments symbol)
  (multiple-value-bind (p found?)
      (gethash (cons symbol arguments) (p d))
    (unless found?
      (warn "Categorical probability of ~a given ~a not found." symbol arguments))
    (pr:in p)))

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
  (let* ((probabilities (mapcar (lambda (s) (probability d arguments s))
				congruent-states))
	 (sum (apply #'pr:add probabilities)))
    (loop for s in congruent-states for p in probabilities do
	 (setf (gethash s table) (pr:div p sum)))
    table))
