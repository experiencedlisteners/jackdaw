(cl:defpackage #:jackdaw
  (:use #:common-lisp)
  (:export
   "GENERATIVE-MODEL" "V" "RECURSIVE" "ACCUMULATOR" "NGRAM-ACCUMULATOR"
   "ONE-SHOT"
   "DEFMODEL" "TRANSITION" "GENERATE-STATES"
   "DISTRIBUTION" "BERNOUILLI" "CATEGORICAL" "UNIFORM"
   "OBSERVE" "PROBABILITY" "PROBABILITIES"
   "+INACIVE+" "+SINGLETON+")
  (:documentation "A toolkit for defining dynamic Bayesian networks 
with congruency constraints."))

(cl:in-package #:jackdaw)

;; Special symbols

(defvar +inactive+ '✳)

;; Globals for enabling meta-programming with models

(defparameter *models* nil) ; list of defined models
(defparameter *model-parameters* nil) ; plist of parameters per model

;; Settings

(defparameter *generate-a-priori-states* nil
  "When T, all a-priori congruent states are generated.
This is useful when we want to calculate things like the entropy of
the predictive distribution.")
(defparameter *estimate?* nil
  "Use to modify the behavior of GENERATE-MOMENT. When T,
GENERATE-MOMENT will only generate values consistent with observations
and will not generate probability distributions")

(defun model-exists? (model-symbol-or-name &optional (package :jackdaw))
  (not (null (find-model model-symbol-or-name package))))

(defun find-model (model-symbol-or-name &optional (package :jackdaw))
  (let ((symbol (find-symbol (if (stringp model-symbol-or-name)
				 (string-upcase model-symbol-or-name)
				 (symbol-name model-symbol-or-name))
			     package)))
    (when (subtypep symbol 'bayesian-network)
      symbol)))

;; Global parameters for CSV output writing.

(defparameter *sequence* nil)
(defparameter *moment* nil)

;; Object serialization helper macros.

(defmacro defreader (type (model data) &body body)
  `(defmethod deserialize ((,model ,type) stream)
     (let ((,data (read stream)))
       ,@body)))

(defmacro defwriter (type (model) &body data)
  `(defmethod serialize ((,model ,type) &optional stream)
     (let ((data (progn ,@data)))
       (unless (null stream) (write data :stream stream))
       data)))

;; Utility

(defmacro %muffle-redefinition-warnings (&body body)
  `(locally
       (declare #+sbcl(sb-ext:muffle-conditions sb-kernel:redefinition-warning))
     (handler-bind
	 (#+sbcl(sb-kernel:redefinition-warning #'muffle-warning))
       (progn ,@body)
      )))

(defun %lambda-list->direct-slots (lambda-list)
  (let ((required)
	(keys)
	(state 'positional)
	(direct-slots))
    (dolist (param lambda-list)
      (case param
	(&key (if (eq state 'positional)
		  (setf state 'key)))
	(t (if (eq state 'positional)
	       (push param required)
	       (if (listp param)
		   (push (list (car param) (cadr param)) keys)
		   (push (list param nil) keys))))))
    (dolist (p required)
      (push `(,p :initarg ,(%kw p) :reader ,p)
		 ;;:initform (%required-arg ,p ,cls))
	    direct-slots))
    (dolist (pdef keys direct-slots)
      (let ((p (car pdef))
	    (default (cadr pdef)))
	(push `(,p :initarg ,(%kw p) :reader ,p :initform ,default) direct-slots)))
    direct-slots))

(defun %lambda-list->params (lambda-list)
  (loop for p in lambda-list
	if (not (member p lambda-list-keywords))
	  collect (if (listp p) (car p) p)))

(defun %lambda-list->plist (lambda-list)
  (apply #'append 
	 (loop for p in (%lambda-list->params lambda-list)
	       collect `(,(%kw p) ,p))))

(defun copy-hash-table (hash-table)
  (let ((ht (make-hash-table 
             :test (hash-table-test hash-table)
             :rehash-size (hash-table-rehash-size hash-table)
             :rehash-threshold (hash-table-rehash-threshold hash-table)
             :size (hash-table-size hash-table))))
    (loop for key being each hash-key of hash-table
       using (hash-value value)
       do (setf (gethash key ht) value)
       finally (return ht))))

(defun %kw (symbol)
  (intern (symbol-name symbol) 'keyword))

(defmacro gethash-or (key hash-table &optional (default-form `(error "~a not found in hash table." ,key)))
  "Like gethash, but default-form is not executed if the key is found.
By default, the default-form will throw an error if the key is not found."
  `(multiple-value-bind (value found?)
       (gethash ,key ,hash-table)
     (if found? value
	 ,default-form)))

;; Objects

(defclass bayesian-network (distribution dag)
  ((%var-specs :allocation :class :type list)
   (%parameter-slots :reader %parameter-slots :allocation :class :type list)
   (output :accessor output :initarg :output :initform nil)
   (output-vars :accessor output-vars :initform nil :initarg :output-vars)
   (distributions :reader distributions :type hash-table)
   (variables :reader variables :type hash-table)))

(defclass dynamic-bayesian-network (bayesian-network) ())

(defclass random-variable ()
  ((name :initarg :name :reader name)
   (output :initarg :output :reader output)
   (distribution :initarg :distribution :accessor distribution :type distribution)
   (model :initarg :model :reader model :type bayesian-network)
   (hidden :initarg :hidden :accessor hidden :initform t)
   (key :initarg :key :reader key)))

(defclass dag ()
  ((vertices :accessor vertices)
   (edge-table :reader edge-table :type 'hastable)))

(defmethod edges ((m dag) vertex)
  (gethash vertex (edge-table m)))

(define-condition dag-contains-cycle (error) ())

(defmethod topological-sort ((graph dag)
			     &optional
			       (vertices (vertices graph))
			       (visited (make-hash-table))
			       result)
  "Sort the vertices of the graph represented by EDGES topologically. Return a list
of vertex indices."
  (if (null vertices) result
      (let ((vertex (car vertices))
	    (remaining (cdr vertices)))
	(topological-sort graph remaining visited (visit graph vertex visited result)))))

(defmethod visit ((graph dag) vertex visited result)
  "A node is a CONS whose CDR is its children and whose CAR is another
CONS the CAR of which is its mark and whose CDR is the a feature index."
  (let ((children (edges graph vertex)))
    (case (gethash vertex visited :no)
      (:yes result)
      (:in-progress
       (error 'dag-contains-cycle))
      (:no
       (progn
	 (setf (gethash vertex visited) :in-progress)
	 (let ((result (topological-sort graph children visited result)))
	   (setf (gethash vertex visited) :yes)
	   (cons vertex result)))))))

(defun hash-table-values (hash-table)
  (loop for value being the hash-values of hash-table collect value))

;; Representation of states and variable identifiers

(defun constr-arg (v)
    "Given a keyword representation of a variable, return
its congruency constraint function symbol. For example,
given :X, return $X"
  (intern (format nil "$~A" (symbol-name v))))

(defun previous (v)
  "Return a symbol referring to the previous value of V in states."
  (intern (format nil "^~A" (symbol-name v)) (symbol-package v)))

(defun horizontal? (v)
  "Return true if V is a horizontal dependency."
  (eq (elt (symbol-name v) 0) #\^))

(defun basename (s)
  "Given an a priori version of a variable name, return its
stem. For example, if S is ^X, (BASENAME S) is X."
  (assert (horizontal? s) ()
	  "Basename ~a not defined, must be horizontal dependency." s)
  (intern (subseq (symbol-name s) 1) (symbol-package s)))

(defmethod horizontal-edges ((m dynamic-bayesian-network) vertex)
  "Return those of dependencies of VERTEX in M that are 
horizontal dependencies."
  (loop for v in (edges m vertex)
     if (horizontal? v)
       collect (basename v)))

(defmethod vertical-edges ((m dynamic-bayesian-network) vertex)
  "Return those dependencies of VERTEX in M that are 
not horizontal."
  (loop for v in (edges m vertex)
     if (horizontal? v)
       collect (basename v)))

(defun get-vertical-arguments (vertices)
  "Return those vertices in VERTICES that are
horizontal dependencies."
  (loop for v in vertices
     if (not (horizontal? v))
     collect v))

(defun inactive? (s)
  "Return T if the value of S is +INACTIVE+."
  (eq s +inactive+))

;; Model definition macro

(defmacro %required-arg (symbol cls)
  (format nil "Slot ~a of ~a should have a value" symbol (type-of cls)))
  ;;`(error ,(format nil "~a is a required initialization argument of ~a."
;;		   (%kw symbol) cls)))

(defun %param-name (p)
  (if (listp p) (car p) p))
  
(defmacro defmodel (class superclasses parameters variables)
  "Model definition macro

Create a model class with name CLASS and superclasses SUPERCLASSES.
PARAMETERS is a lambda list which supports default values and &key 
arguments but not &optional or (a default a-supplied-p) style parameters.
VARIABLES is a list of variable definitions."
  (let* ((direct-slots (%lambda-list->direct-slots parameters))
	 (parameter-names (mapcar #'%param-name (remove '&key parameters)))
	 (edges (make-hash-table))
	 (dist-specs) (var-specs) (vertices)
	 (methods))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c) (subtypep c 'jackdaw:bayesian-network)) superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type BAYESIAN-NETWORK." class)
    ;; Initialize dist-specs, var-specs, and vertices, verify model consistency,
    ;; and define congruency constraint methods.
    (loop
      for variable in (reverse variables) ; since we're PUSHing
      collect
      (destructuring-bind (v parents distribution constraint
			   &key (key `(lambda (moment)
					(getf moment ,(%kw v)
					      (format nil "~s not found in moment" ',v))))
			     (formatter '#'identity) (hidden t))
	  variable
	(setf (gethash v edges) parents)
	(push v vertices)
	(push distribution dist-specs)
	(push (list key formatter hidden) var-specs)
	(push 
	 `(defmethod congruent-values ((model ,class)
				      (variable (eql ',v))
				      parent-state)
	    (declare (ignorable model))
	    (let* (,@(loop for p in parents
			   collect `(,(constr-arg p) (gethash ',p parent-state)))
		   ,@(loop for parameter in parameter-names collect
			   `(,parameter (slot-value model ',parameter)))
		   ,@(when (member (previous v) parents)
		       `(($^self ,(constr-arg (previous v))))))
	      (declare (ignorable
			,@(mapcar #'constr-arg parents)
			,@parameter-names
			,@(when (member (previous v) parents) `($^self))))
	      (handler-case ,constraint
		(error (e)
		  (warn
		   "An error occurred in the a priori congruency constraint of ~A." ',v)
		  (error e)))))
	 methods)))
    `(%muffle-redefinition-warnings
       (pushnew (%kw ',class) *models*)
       (setf (getf *model-parameters* ',class) ',parameters)
       (defclass ,class ,(if (null superclasses) '(bayesian-network) superclasses)
	 ((%var-specs :initform ',var-specs)
	  (%parameter-slots :initform ',parameter-names)
	  ;;(%dist-specs :initform ',dist-specs)
	  ,@direct-slots))
       ,@methods
       (defmethod vertices ((m ,class))
	 ',vertices)
       (defmethod edge-table ((m ,class)) ,edges)
       (defmethod initialize-instance :after ((model ,class) &key)
	 (let (,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
			collect `(,p (,p model))))
	   (declare (ignorable ,@parameter-names))
	    ,@(loop for v in vertices
		    for dist-spec in dist-specs
		    collect
		    (destructuring-bind (dist dist-args &rest dist-params) dist-spec
		      (assert (subsetp dist-args (gethash v edges))
			      () "Arguments to distribution of ~a must be subset of parents." v)
		      (assert (subtypep dist 'distribution) ()
			      "Distribution argument of the variable ~a (~a) must be of type DISTRIBUTION." 
			      v dist)
		      `(setf (distribution (gethash ',v (variables model)))
			     (apply #'make-instance ',dist
				    ,(append `(list :arguments ',dist-args
						    :variable-symbol ',v)
					     dist-params)))))))
       (defun ,(intern (format nil "MAKE-~A-MODEL" (symbol-name class)))
	   (,@parameters ,@(unless (member '&key parameters) '(&key)) output output-vars observe)
	 (make-instance ',class :output output :output-vars output-vars :observe observe
			,@(%lambda-list->plist parameters))))))

(defmethod initialize-instance :after ((model bayesian-network) &key observe)
  (dolist (v (output-vars model))
    (assert (member v (vertices model)) ()
	    "Output variable ~a is not a variable of ~a." v (type-of model)))
  (let ((variables (make-hash-table)))
    (loop for v in (vertices model)
	  for (key formatter hidden) in (slot-value model '%var-specs)
	  ;;for dist-spec in (slot-value model '%dist-specs)
	  do
	     (setf (gethash v variables)
		   (make-instance
		    'jackdaw::random-variable
		    :model model
		    :hidden hidden
		    :name v :key (eval key) :output (eval formatter))))
    (setf (slot-value model 'variables) variables))
  (apply #'observe (cons model observe)))

;; Model serialization

(defwriter random-variable (v)
  (serialize (distribution v)))
(defreader random-variable (v data)
  (deserialize (distribution v) data))
(defwriter bayesian-network (m)
  (let ((variables)
	(parameter-slots (loop for s in (%parameter-slots m)
			       collect (slot-value m s))))
    (loop for vertex being the hash-key using (hash-value variable) of (variables m) do
      (setf (getf variables variables) (serialize variable)))
    (list parameter-slots variables)))
(defreader bayesian-network (m data)
  (destructuring-bind (parameters variables)
      data
    (loop for p in parameters
	  for s in (%parameter-slots m)
	  do (setf (slot-value m s) p))
    (loop for (vertex variable) on variables by #'cddr do
      (with-input-from-string (s (write-to-string variable))
	(deserialize (gethash vertex (variables m)) s)))))

;; Jackdaw mechanics

(defmethod formatter-function ((m bayesian-network) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((v (gethash variable (variables m))))
    (output v)))

(defmethod observed-value ((m bayesian-network) variable moment)
  "Obtain the observed value of VARIABLE from MOMENT given M."
  (let ((v (gethash variable (variables m))))
    (handler-case
	(funcall (key v) moment)
      (error (e)
	(warn "An error occurred in the observation function of ~a"
	      variable)
	(error e)))))

(defmethod hidden? ((m bayesian-network) variable)
  (hidden (gethash variable (variables m))))

(defmethod observed? ((m bayesian-network) variable)
  (not (hidden (gethash variable (variables m)))))

(defmethod %set-hidden ((m bayesian-network) hidden vertices)
  (dolist (v vertices)
    (let ((variable
	    (gethash-or v (variables m)
		       (error "~a is not a variable of ~a" v (type-of m)))))
      (setf (hidden variable) hidden))))

(defmethod hide ((m bayesian-network) &rest vertices)
  "Hide VARIABLES."
  (%set-hidden m t vertices))
  
(defmethod observe ((m bayesian-network) &rest vertices)
  "Make VARIABLES observable."
  (%set-hidden m nil vertices))

(defmethod model-variable-distribution ((m bayesian-network) vertex)
  (distribution (gethash vertex (variables m))))

(defmethod state-variables ((m dynamic-bayesian-network))
  (reduce #'union (mapcar (lambda (v) (horizontal-edges m v)) (vertices m))))

(defmethod model-variables ((m dynamic-bayesian-network))
  (union (state-variables m) (observed-variables m)))

(defun get-probability (state distribution)
  (let ((probability (gethash-or state distribution
			 (warn "State ~a not found in distribution." state))))
    probability))

(defmethod order ((m bayesian-network))
  (length (vertices m)))

(defmethod observations ((m bayesian-network) moment)
  (let ((observations (make-hash-table)))
    (dolist (v (observed-variables m) observations)
      (setf (gethash v observations) (observed-value m v moment)))))

(defmethod congruent-states (states observations)
  "Return the subset of STATES that is consistent with OBSERVATIONS."
  (loop for s in states if (observed? s observations) collect s))

(defmethod congruent? (state observations)
  "Return T if STATE is consistent with OBSERVATIONS."
  (every #'identity
	 (loop for v being each hash-key of observations
	       collect (equal (gethash v observations) (gethash-or v state)))))

(defmethod observed-variables ((m bayesian-network))
  (remove-if (lambda (v) (hidden? m v)) (vertices m)))

(defmethod make-root-state ((m bayesian-network))
  "Every root node in a Bayesian network is implicitly conditioned
on this root state."
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) (pr:in 1))
    state))

(defmethod make-root-state ((m dynamic-bayesian-network))
  "Every root node in a Bayesian network is implicitly conditioned
on this root state."
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) (pr:in 1))
    (dolist (variable (state-variables m))
      (setf (gethash variable state) +inactive+))
    state))

(defun branch-state (probability variable value &optional origin)
  (let ((new-state (if (null origin) (make-hash-table)
		       (copy-hash-table origin))))
    (setf (gethash :probability new-state) probability)
    (setf (gethash variable new-state) value)
    new-state))


(defun format-obs (observation)
  (format nil "(~{~{~w~^: ~}~^, ~})~%"
	  (loop for k being each hash-key of observation
		collect (list k (gethash k observation)))))

(defmethod generate-observations ((m generative-model) dataset)
  "Generate a PLIST with for each variable a list containing the 
variable's value in each moment in each sequence. For this to work,
M must be fully observed in the sense that in any moment,
each variable has exactly one possible value."
  (let* ((*estimate?* t)
	 (*generate-a-priori-states* nil)
	 (observation-dataset))
    (dolist (sequence dataset)
      (let ((state (make-root-state m))
	    (observation-sequence))
	(dolist (moment sequence)
	  (let* ((observation (observations m moment))
		 (states (generate-moment m observation
					  (rotate-state m state :keep-trace? nil))))
	    (assert (> (length states) 0) ()
		    "Could not generate observation ~a of moment ~w" (format-obs observation) moment)
	    (assert (eq (length states) 1) ()
		    "Cannot estimate from model that is not fully observed.")
	    	    (push (loop for v in (vertices m) collect (gethash v (car states)))
		  observation-sequence)
	    (setf state (car states))))
	(push (reverse observation-sequence) observation-dataset)))
    (reverse observation-dataset)))

(defmethod estimate ((m dynamic-bayesian-network) dataset)
  (warn "generating observations")
  (let* ((dataset (generate-values m dataset)))
    (dotimes (vertex-index (order m) m)
      (let* ((v (elt (vertices m) vertex-index))
	     (distribution (model-variable-distribution m v))
	     (root (loop repeat (order m) collect +inactive+))
	     (arguments
	       (mapcar (lambda (v)
			 (position (if (horizontal? v) (basename v) v) (vertices m)))
		       (arguments distribution)))
	     (horizontal?
	      (mapcar (lambda (v) (horizontal? v)) (arguments distribution))))
	(labels
	    ((get-variable-observation (previous current)
	       (mapcar (lambda (i horizontal?) (elt (if horizontal? previous current) i))
		       (cons vertex-index arguments)
		       (cons nil horizontal?)))
	     (get-observation-sequence (sequence)
	       (loop for (previous current) on (cons root sequence) by #'cdr
		     repeat (length sequence)
		     collect (get-variable-observation previous current))))
	  (warn "converting data to observations of ~a" v)
	  (let ((dataset (mapcar #'get-observation-sequence dataset)))
	    (warn "estimating ~a" v)
	    ;;(print v)
	    ;;(print dataset)
	    (estimate distribution dataset)))))))

(defmethod probability ((m bayesian-network) observation)
  (evidence m (generate m observation)))

(defmethod probabilities ((d distribution) parents-state congruent-values)
  "Obtain the probability of provided CONGRUENT-VALUES by the PROBABILITY method.
This is just a wrapper for PROBABILITY-DISTRIBUTION which grabs arguments from PARENTS-STATE
and avoids a call to PROBABILITY-DISTRIBUTION when the variable is inactive."
  (let ((arguments (mapcar (lambda (v) (gethash v parents-state)) (arguments d))))
    (if (equal congruent-values (list +inactive+))
	(list (pr:in 1))
	(probability-distribution d arguments congruent-values))))

(defmethod probability-distribution ((d distribution) arguments congruent-values)
  (mapcar (lambda (s) (probability d (cons s arguments)))
	  congruent-values))

(defmethod descr ((m bayesian-network))
  (format nil "~a model" (symbol-name (type-of m))))

(defmethod descr ((v random-variable))
  (format nil "variable ~w of ~a" (name v) (descr (model v))))

(defmethod descr ((d distribution))
  (format nil "~a distribution of variable ~a" (type-of d) (variable-symbol d)))

;;(defmethod generate :before (v obs &optional p)
;;  (format t "Generating ~a~%" (descr v)))

(defmethod generate ((variable random-variable) observation
		     &optional (parent-state (make-hash-table)))
  (let* ((vertex (name variable))
	 (distribution (distribution variable))
	 (a-priori-congruent
	   (congruent-values (model variable) vertex parent-state))
	 (hidden? (hidden variable))
	 (parent-probability (gethash :probability parent-state))
	 (congruent-values
	   (if (and (not hidden?))
	       (when (member (gethash-or vertex observation)
			     a-priori-congruent)
		 (list (gethash vertex observation)))
	       a-priori-congruent)))
    (flet ((normalize (probabilities)
	     (let ((total-mass (apply #'pr:add probabilities)))
	       (mapcar (lambda (p) (pr:div p total-mass)) probabilities)))
	   (branch-state (value &optional probability)
	     (branch-state (unless (null probability)
			     (pr:mul probability parent-probability))
			   vertex value parent-state)))
      (if *estimate?*
	  (mapcar #'branch-state congruent-values)
	  (let* ((probabilities (probabilities distribution parent-state congruent-values))
		 (probabilities (if hidden? (normalize probabilities)
				    probabilities)))
	    (mapcar #'branch-state congruent-values probabilities))))))

(defmethod generate ((m bayesian-network) observation
		     &optional (parent-state (make-root-state m)))
  (generate-moment m parent-state observation))

(defmethod generate ((m dynamic-bayesian-network) observation
		     &optional (parent-state (list (make-root-state m))))
  (generate-sequence m observation :congruent-states parent-state))
		     
(defun generate-moment (m observation &optional (parent-state (make-root-state m))
					(vertices (get-vertical-arguments (topological-sort m))))
  "Generate the congruent states that can be reached from PARENT-STATE."
  (let* ((parent-states
	   (if (null (cdr vertices)) (list parent-state)
	       (generate-moment m observation parent-state (cdr vertices))))
	 (new-states))
    (dolist (parent-state parent-states new-states)
      (let ((variable (gethash (car vertices) (variables m))))
	(dolist (state (generate variable observation parent-state))
	  (push state new-states))))))

(defun make-marginal-state (variables values trace)
  (let ((state (make-hash-table)))
    (setf (gethash :trace state) trace)
    (loop for variable in variables for value in values do
      (setf (gethash variable state) value))
    state))

(defun update-marginal-state (marginal-state p)
  (let ((state-p (gethash :probability marginal-state)))
    (setf (gethash :probability marginal-state)
	  (apply #'pr:add (cons p (when state-p (list state-p)))))
    marginal-state))

(defun marginalize (states variables)
  "Marginalize a list of states wrt. variables and return a new list of states."
  (let ((marginal (make-hash-table :test #'equal)))
    (dolist (state states)
      (let* ((trace (gethash :trace state))
	     (values (loop for v in variables collect (gethash v state)))
	     (marginal-state (gethash-or values marginal
				 (make-marginal-state variables values
						      trace)))
	     (probability (gethash :probability state))
	     (new-state (update-marginal-state marginal-state probability)))
	(setf (gethash values marginal) new-state)))
    (hash-table-values marginal)))

(defmethod rotate-state ((m dynamic-bayesian-network) state &key (keep-trace? t))
  "\"Rotate\" a state. In  rotated (a priori) version of a state, every parameter
X is renamed ^X and variables of the form ^X in STATE are dropped.
Optionally, a trace can be kept which stores a hash table containing the dropped
values, their probability and the previous trace."
  (let ((new-state (make-hash-table))
	(persistent-variables (model-variables m)))
    (setf (gethash :probability new-state) (gethash :probability state))
    (when keep-trace?
      (let ((trace (make-hash-table)))
	(dolist (key (cons :trace persistent-variables))
	  (setf (gethash key trace) (gethash key state)))
	(setf (gethash :trace new-state) trace)))
    (dolist (variable persistent-variables new-state)
      (setf (gethash (previous variable) new-state)
	    (gethash variable state)))))
  
(defmethod transition ((m dynamic-bayesian-network) moment congruent-states
		       &key keep (keep-trace? t))
  (let ((observation (observations m moment))
	(persistent-variables (union (model-variables m) keep)))
    (flet ((generate-moment (state)
	     (let* ((states
		      (generate-moment m observation
				       (rotate-state m state :keep-trace? keep-trace?)))
		    (congruent-states
		      (if *generate-a-priori-states*
			  (congruent-states states observation)
			  states)))
	       (%write-states m congruent-states observation)
	       (marginalize congruent-states persistent-variables))))
      (let* ((new-states
	       (apply #'append (mapcar #'generate-moment congruent-states))))	
	(when (eq (length new-states) 0)
	  (error "No a posteriori congruent states at moment
~a at position #~a in sequence ~a." moment *moment* *sequence*))
	(marginalize new-states persistent-variables)))))

(defmethod generate-sequences ((m dynamic-bayesian-network) dataset &optional (write-header? t))
  "Process a dataset of sequences. Note that this method assumes each sequence is a CONS 
the CAR of which is used to identify the sequence. This may be something like an index 
or a unique-id."
  (let* ((sequence (car dataset))
	 (*sequence* (car sequence)))
    ;;(warn "~a" (car sequence))
    (generate-sequence m (cdr sequence) :write-header? write-header?))
  (unless (null (cdr dataset))
    (generate-sequences m (cdr dataset) nil)))

(defmethod generate-sequence ((m dynamic-bayesian-network) moments
			  &key keep
			    (keep-trace? t)
			    (write-header? t)
			    (moment 0)
			    (congruent-states (list (make-root-state m))))
  "Process a sequence moment by moment and return the states that remained 
congruent by the end of the sequence."
  (when write-header? (%write-header m))
  (if (null moments)
      (dolist (v (vertices m) congruent-states)
	(next-sequence (model-variable-distribution m v) congruent-states))
      (let* ((*moment* moment)
	     (new-congruent-states (transition m (car moments) congruent-states
					       :keep keep :keep-trace? keep-trace?)))
	;;(format t "Evidence ~a, prob first con state: ~a n-cong: ~a~%"
	;;	(evidence m new-congruent-states)
	;;	(gethash :probability (car new-congruent-states))
	;;	(length new-congruent-states))
	(generate-sequence m (cdr moments)
			   :keep keep
			   :write-header? nil
			   :moment (1+ moment)
			   :congruent-states new-congruent-states))))

(defmethod evidence ((m bayesian-network) congruent-states)
  (gethash :probability
	   (car (marginalize congruent-states (observed-variables m)))))

(defmethod posterior-distribution ((m bayesian-network) congruent-states)
  (let ((evidence (evidence m congruent-states)))
    (dolist (state congruent-states congruent-states)
      (setf (gethash :probability state) (pr:div (gethash :probability state) evidence)))))

(defun states->probabilities (states &rest variables)
  (loop for s in states
	collect
	(list (loop for v in variables collect (gethash v s))
	      (gethash :probability s))))
			      
(defmethod %write-header ((m bayesian-network) &optional (output (output m)))
  (let* ((output-var-names
	   (loop for v in (if (null (output-vars m)) (vertices m)
			      (output-vars m))
		 collect (string-downcase (symbol-name v))))
	 (columns (append
		   '("sequence" "moment")
		   output-var-names
		   (when *generate-a-priori-states* '("congruent"))
		   '("probability"))))
    (format output "~{~a~^,~}~%" columns)))

(defmethod %write-states ((m bayesian-network) states observations &optional (output (output m)))
  (dolist (state states)
    (%write-state m state observations output)))

(defmethod %write-state ((m bayesian-network) state observations &optional (output (output m)))
  (format output "~a,~a~{,~a~^~},~a,~a~%" *sequence* *moment*
	  (loop for v in (if (null (output-vars m))
			     (vertices m)
			     (output-vars m))
		collect (funcall (formatter-function m v) (gethash v state)))
	  (congruent? state observations) (gethash :probability state)))

(defun trace-back (state variable &optional trace)
  (let ((new-trace (cons (gethash variable state) trace))
	(previous-state (gethash :trace state)))
    (if (null previous-state)
	new-trace
	(trace-back previous-state variable new-trace))))

(defmethod pprint-state ((m bayesian-network) state)
  (format t "(~{~{~A~^:~}~^ ~}): ~a"
	  (loop for k being each hash-key of state
		if (not (eq k :probability))
		  collect (list k (gethash k state)))
	  (gethash :probability state)))
