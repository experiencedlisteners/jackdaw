(cl:defpackage :jackdaw
  (:use #:common-lisp)
  (:nicknames :jd)
  (:export
   ;; Probability distributions
   #:probability-distribution
   #:conditional-probability-distribution
   #:ngram-model #:make-ngram-model-distribution
   #:bernouilli #:make-bernouilli-distribution
   #:cpt #:make-cpt-distribution
   #:ppms #:make-ppms-distribution
   #:uniform #:make-uniform-distribution
   ;; Probability distributions API
   #:probability #:conditional-probability
   #:conditional-probabilities #:estimate
   ;; Bayesian networks
   #:observation
   #:bayesian-network #:dynamic-bayesian-network
   #:random-variable
   ;; Evaluating models
   #:generate
   #:generate-sequence ;; TODO REMOVE
   #:sequence-generate #:probability
   #:generate-congruent-values
   #:probabilities #:domain #:estimate
   ;; Model properties
   #:variables #:observed-variables #:model-variable
   #:model-variable-distribution
   ;; Variable properties
   #:vertex #:distribution #:parents
   ;; Manipulating Bayesian networks
   #:observe #:hide
   ;; Manipulating lists of states
   #:transition #:evidence #:posterior
   #:state-probability-table #:trace-back
   #:marginalize
   #:pprint-state 
   ;; Model definition tools
   #:defmodel #:defdistribution #:defestimator
   ;; Congruency constraint utilities
   #:deterministic #:normal
   #:recursive #:persist #:one-shot
   #:accumulate #:ngram #:chain #:markov
   #:ensure-list #:inactive?
   ;; Serialization utilities
   #:defreader #:defwriter
   ;; Metaprogramming
   #:model-exists? #:find-model
   ;; Constants
   #:+inactive+
   ;; Settings
   #:*models* #:*estimate?* #:*model-parameters*
   #:*generate-a-priori-states*)
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
      (push `(,p :initarg ,(%kw p) :accessor ,p)
		 ;;:initform (%required-arg ,p ,cls))
	    direct-slots))
    (dolist (pdef keys direct-slots)
      (let ((p (car pdef))
	    (default (cadr pdef)))
	(push `(,p :initarg ,(%kw p) :accessor ,p :initform ,default) direct-slots)))
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

(defclass dag ()
  ((vertices :accessor vertices :initform nil)
   (edge-table :reader edge-table :type 'hastable)))

(defclass probability-distribution ()
  ((%parameters :reader %parameters :initform nil)))

(defclass conditional-probability-distribution (probability-distribution) ())

(defclass bayesian-network (probability-distribution dag)
  ((%var-specs :allocation :class :type list :initform nil)
   (output :accessor output :initarg :output :initform nil)
   (output-vars :accessor output-vars :initform nil :initarg :output-vars)
   (distributions :reader distributions :type hash-table)
   (variables :reader variables :type hash-table)))

(defclass dynamic-bayesian-network (bayesian-network)
  ((hidden-state :accessor hidden-state)
   (keep-trace? :initarg :keep-trace? :accessor keep-trace? :initform nil)
   (intermediate-marginalization :initarg :intermediate-marginalization
				 :accessor intermediate-marginalization)))

(defclass random-variable ()
  ((parents :initarg :parents :reader parents :initform nil)
   (vertex :initarg :vertex :reader vertex)
   (output :initarg :output :reader output)
   (distribution :initarg :distribution :accessor distribution :type probability-distribution)
   (distribution-parents :initarg :distribution-parents :accessor distribution-parents)
   (model :initarg :model :reader model :type bayesian-network)
   (hidden :initarg :hidden :accessor hidden :initform t)
   (observer :initarg :observer :accessor observer)))

(defclass previous-random-variable ()
  ((vertex :initarg :vertex :reader vertex)
   (model :initarg :model :reader model :type dynamic-bayesian-network)))

(defmethod edges ((m dag) vertex)
  (gethash vertex (edge-table m)))

(defmethod vertices ((m dynamic-bayesian-network))
  (append (mapcar #'previous (state-variables m))
	  (slot-value m 'vertices)))

(defmethod order ((m dag))
  (length (vertices m)))

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
  (if (not (horizontal? s)) s
      (intern (subseq (symbol-name s) 1) (symbol-package s))))

(defmethod horizontal-edges ((m dynamic-bayesian-network) vertex)
  "Return those of dependencies of VERTEX in M that are 
horizontal dependencies."
  (loop for v in (edges m vertex)
     if (horizontal? v)
       collect v))

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
	 (vertices (mapcar #'car variables))
	 (observer-arguments
	   (loop for v in vertices collect (intern (format nil "~a-OBSERVER" v))))
	 (edges (make-hash-table))
	 (dist-specs) (var-specs) (methods))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c)
					 (subtypep c 'jackdaw:bayesian-network))
				       superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type BAYESIAN-NETWORK." class)
    ;; Initialize dist-specs, var-specs, and vertices, verify model consistency,
    ;; and define congruency constraint methods.
    (loop
      for variable in (reverse variables) ; since we're PUSHing
      collect
      (destructuring-bind (v parents distribution constraint
			   &key (observer `(lambda (moment)
					(getf moment ,(%kw v)
					      (format nil "~s not found in moment" ',v))))
			     (formatter '#'identity) (hidden t))
	  variable
	(setf (gethash v edges) parents)
	(assert (subsetp (mapcar #'basename parents) vertices) ()
		"The parents {~{~a~^, ~}} of ~a are not known variables of ~a"
		(set-difference (mapcar #'basename parents) vertices)
		v class)
	
	(push distribution dist-specs)
	(push (list observer formatter hidden) var-specs)
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
	  (%parameters :initform ',parameter-names)
	  (vertices :initform ',vertices)
	  ;;(%dist-specs :initform ',dist-specs)
	  ,@direct-slots))
       ,@methods
       (defmethod edge-table ((m ,class)) ,edges)
       (defmethod initialize-instance :after
	   ((model ,class)
	    &key ,@observer-arguments)
	 (let (,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
			collect `(,p (,p model))))
	   (declare (ignorable ,@parameter-names))
	   ,@(loop for v in vertices for arg in observer-arguments 
		   collect
		   `(unless (null ,arg)
		      (setf (observer (model-variable model ',v)) ,arg)))
	   ,@(loop for v in vertices
		   for dist-spec in dist-specs
		   collect
		   (destructuring-bind (dist dist-args &rest dist-params) dist-spec
		     (assert (subsetp dist-args (gethash v edges))
			     () "Arguments to distribution of ~a must be subset of parents." v)
		     (assert (subtypep dist 'probability-distribution) ()
			     "Distribution argument of the variable ~a (~a) must be of type PROBABILITY-DISTRIBUTION." 
			     v dist)
		     `(setf (distribution-parents (model-variable model ',v))
			    ',dist-args
			    (distribution (model-variable model ',v))
			    (apply #'make-instance ',dist (list ,@dist-params)))))))
       (defun ,(intern (format nil "MAKE-~A-MODEL" (symbol-name class)) (symbol-package class))
	   (,@parameters ,@(unless (member '&key parameters) '(&key)) output output-vars
	    observe)
	 (make-instance ',class :output output :output-vars output-vars :observe observe
			,@(%lambda-list->plist parameters))))))

(defmethod initialize-instance :before ((model dynamic-bayesian-network) &key observe)
  (dolist (v observe)
    (assert (not (horizontal? v)) ()
	    "Cannot observe ~a in ~a because it is a horizontal dependency."
	    v (type-of model))))

(defmethod initialize-instance :after ((model bayesian-network) &key observe)
  (dolist (v (output-vars model))
    (assert (member v (vertices model)) ()
	    "Output variable ~a is not a variable of ~a." v (type-of model)))
  (let ((variables (make-hash-table)))
    (loop for v in (slot-value model 'vertices)
	  for (observer formatter hidden) in (slot-value model '%var-specs)
	  do
	     (setf (gethash v variables)
		   (make-instance
		    'jackdaw::random-variable
		    :model model
		    :parents (edges model v)
		    :hidden hidden
		    :vertex v :observer (eval observer) :output (eval formatter))))
    (setf (slot-value model 'variables) variables))
  (unless (null observe) (apply #'observe (cons model observe))))

(defmethod initialize-instance :after ((model dynamic-bayesian-network) &key)
  (let ((variables (variables model)))
    (loop for v in (state-variables model)
	  do
	     (setf (gethash (previous v) variables)
		   (make-instance
		    'previous-random-variable
		    :model model
		    :vertex v)))))

;; Model serialization

(defwriter probability-distribution (m)
  (loop for s in (%parameters m) collect (slot-value m s)))
(defreader probability-distribution (m data)
  (loop for v in data for s in (%parameters m)
	collect (slot-value m s)))
(defwriter random-variable (v)
  (serialize (distribution v)))
(defreader random-variable (v data)
  (deserialize (distribution v) data))
(defwriter bayesian-network (m)
  (let ((variables)
	(parameter-slots (loop for s in (%parameters m)
			       collect (slot-value m s))))
    (loop for vertex being the hash-key using (hash-value variable) of (variables m) do
      (setf (getf variables variables) (serialize variable)))
    (list parameter-slots variables)))
(defreader bayesian-network (m data)
  (destructuring-bind (parameters variables)
      data
    (loop for p in parameters
	  for s in (%parameters m)
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
	(funcall (observer v) moment)
      (error (e)
	(warn "An error occurred in the observation function of ~a"
	      variable)
	(error e)))))


(defmethod observation ((m dynamic-bayesian-network) moment &rest variables)
  (dolist (v variables)
    (assert (not (horizontal? v)) ()
	    "~a represents the past and cannot be observed"
	    v))
  (call-next-method))

(defmethod observation ((m bayesian-network) input &rest variables)
  "Generate a partially observed state of M."
  (dolist (v variables)
    (assert (member v (vertices m)) ()
	    "~a is not a variable of ~a" v (type-of m)))
  (let ((observation (make-hash-table))
	(variables (if (null variables) (observed-variables m)
		       variables)))
    (setf (gethash :probability observation) (pr:in 1))
    (when *generate-a-priori-states*
      (setf (gethash :congruent? observation) t))
    (dolist (v variables observation)
      (setf (gethash v observation) (observed-value m v input)))))

(defmethod hidden? ((m bayesian-network) vertex)
  (hidden (model-variable m vertex)))

(defmethod hidden? ((m dynamic-bayesian-network) vertex)
  (if (horizontal? vertex) t
      (call-next-method)))

(defmethod observed? ((m bayesian-network) variable)
  (not (hidden (gethash variable (variables m)))))

(defmethod %set-hidden ((m bayesian-network) hidden vertices)
  (let ((vertices (if (null vertices) (vertices m) vertices)))
    (dolist (v vertices)
      (let ((variable
	      (gethash-or v (variables m)
		  (error "~a is not a variable of ~a" v (type-of m)))))
	(setf (hidden variable) hidden)))))

(defmethod hide ((m bayesian-network) &rest vertices)
  "Hide VARIABLES."
  (%set-hidden m t vertices)
  (format *error-output*
	  "The following variables of ~a are now observed: ~{~a~^, ~}~%"
	  (type-of m)
	  (observed-variables m)))
  
(defmethod observe ((m bayesian-network) &rest vertices)
  "Make VARIABLES observable."
  (%set-hidden m nil vertices)
  (format *error-output*
	  "The following variables of ~a are now observed: ~{~a~^, ~}~%"
	  (type-of m)
	  (observed-variables m)))

(defmethod model-variable ((m bayesian-network) vertex)
  (gethash-or vertex (variables m)))

(defmethod model-variable-distribution ((m bayesian-network) vertex)
  (distribution (model-variable m vertex)))

(defmethod state-variables ((m dynamic-bayesian-network))
  (remove-duplicates
   (apply #'append
	  (mapcar (lambda (v) (mapcar #'basename (horizontal-edges m v)))
		  (slot-value m 'vertices)))))

(defmethod model-variables ((m dynamic-bayesian-network))
  (union (state-variables m) (observed-variables m)))

(defun get-probability (state distribution)
  (let ((probability (gethash-or state distribution
			 (warn "State ~a not found in distribution." state))))
    probability))

(defmethod congruent-states (states)
  "Return the subset of STATES that has been marked congruent."
  (loop for s in states if (congruent? s) collect s))

(defmethod congruent? (state)
  "Return T if STATE has been marked congruent."
  (if *generate-a-priori-states*
      (gethash-or :congruent? state)
      t))

(defmethod observed-variables ((m bayesian-network))
  (remove-if (lambda (v) (hidden? m v)) (slot-value m 'vertices)))

(defmethod make-root-state ((v random-variable))
  "Every root node in a Bayesian network is implicitly conditioned
on this root state."
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) (pr:in 1))))

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
    (dolist (vertex (state-variables m))
      (setf (gethash vertex state) +inactive+))
    state))

(defun branch-state (probability variable value origin)
  (let ((new-state (copy-hash-table origin)))
    (setf (gethash :probability new-state) probability)
    (setf (gethash variable new-state) value)
    (when *generate-a-priori-states*
      (setf (gethash :congruent? new-state)
	    (equal value (gethash variable origin))))
    new-state))

(defun format-obs (observation)
  (format nil "(~{~{~w~^: ~}~^, ~})~%"
	  (loop for k being each hash-key of observation
		collect (list k (gethash k observation)))))

;; Input, bayesnet -> list of possible outcomes
(defmethod generate-values ((m bayesian-network) input &optional vertices)
  (let ((vertices (or vertices (vertices m)))
        (*estimate?* t)
	(*generate-a-priori-states* nil)
	(states (generate m input)))
    (loop for s in states
	  collect
	  (loop for v in vertices collect (gethash-or v s)))))

;; Input, dyn bayesnet -> list of lists of possible outcomes
;; Technically, the Cartesian-product of these are the possible outcomes of
;; the dynamic Bayesan network.
(defmethod generate-values ((m dynamic-bayesian-network) input-sequence
			 &optional vertices)
  (loop for input in input-sequence
	collect
	(prog1
	    (let ((states (call-next-method m input vertices)))
	      (setf (hidden-state m) (transition-states m states))
	      states))))

(defmacro first-and-only
    (l &optional (message "First and only requires that list contains exactly one item"))
  `(progn
     (assert (eq (length ,l) 1) () ,message)
     (first ,l)))

;; Input, BN -> possible outcome (= outcome of fully observed BN)
(defmethod fully-observed-value ((m bayesian-network) observation &optional vertices)
  (let* ((congruent-values (generate-values m observation vertices)))
    (first-and-only congruent-values
	"Cannot estimate from model that is not fully observed.")))

;; Input, DBN -> list of outcome each moment (= outcome of fully observed DBN)
(defmethod fully-observed-value ((m dynamic-bayesian-network) observation &optional vertices)
  (unless (null observation)
    (cons (fully-observed-value m (car observation) vertices)
	  (fully-observed-value m (cdr observation) vertices))))

;; Inputs, BN -> list of outcomes
(defmethod fully-observed-dataset ((m bayesian-network) observations &optional vertices)
  (loop for observation in observations
	collect	(fully-observed-value m observation vertices)))

(defmethod observed-outcome->variable-outcome ((m bayesian-network) outcome indices)
  (mapcar (lambda (i) (elt outcome i)) indices))

(defmethod observed-outcome->variable-outcome ((m dynamic-bayesian-network) outcome indices)
  (mapcar (lambda (moment) (call-next-method m moment indices)) outcome))
 
(defmethod estimate ((m bayesian-network) dataset)
  (format *error-output* "Generating observations~%")
  (let* ((dataset (fully-observed-dataset m dataset)))
    (dotimes (vertex-index (order m))
      (let* ((vertex (elt (vertices m) vertex-index))
	     (variable (model-variable m vertex))
	     (argument-indices
	       (mapcar (lambda (v) (position v (vertices m)))
		       (distribution-parents variable)))
	     (distribution (model-variable-distribution m vertex)))
	(format *error-output* "Converting data to observations of ~a~%" vertex)
	(let ((dataset (mapcar (lambda (observation)
				 (observed-outcome->variable-outcome
				  m observation argument-indices))
			       dataset)))
	  (format *error-output* "Estimating ~a~%" vertex)
	  (estimate distribution dataset))))))

(defmethod probability ((m bayesian-network) input)
  (evidence m (generate m (observation m input))))

(defmethod probability ((m dynamic-bayesian-network) input)
  (sequence-generate m input :marginal (observed-variables m))
  (evidence m (hidden-state m)))

(defmethod conditional-probabilities ((variable random-variable) congruent-values
				      &optional parents-state)
  "Obtain the probabilities of a list of CONGRUENT-VALUES of VARIABLE.
This is just a wrapper for the PROBABILITIES of the variable's distribution.
It grabs the arguments from the parents
and avoids a call to PROBABILITY-DISTRIBUTION when the variable is inactive."
  (let ((arguments (mapcar (lambda (v) (gethash v parents-state)) (distribution-parents variable))))
    (if (equal congruent-values (list +inactive+))
	(list (pr:in 1))
	(conditional-probabilities (distribution variable)
				   congruent-values arguments))))

(defmethod descr ((m bayesian-network))
  (format nil "~a model" (symbol-name (type-of m))))

(defmethod descr ((v random-variable))
  (format nil "variable ~w of ~a" (vertex v) (descr (model v))))

(defmethod descr ((d probability-distribution))
  (format nil "~a" (type-of d)))

(defmethod initialize ((m dynamic-bayesian-network))
  (setf (hidden-state m) (list (make-root-state m))))

;; model, state -> list of states
(defmethod graph-generate ((m bayesian-network) parent-state
			   &optional (sorted-vertices (topological-sort m)))
  (let* ((parent-states
	   (if (null (cdr sorted-vertices)) (list parent-state)
	       (graph-generate m parent-state (cdr sorted-vertices))))
	 (new-states))
    (dolist (parent-state parent-states new-states)
      (let ((variable (model-variable m (car sorted-vertices))))
	(dolist (state (generate variable parent-state))
	  (push state new-states))))))

;;(defmethod generate :before (v obs &optional p)
;;  (format t "Generating ~a~%" (descr v)))

(defmethod generate ((variable previous-random-variable)
		     &optional (parent-state (make-root-state variable)))
  "Variables in the previous moment are already observed and their probability
is already incorporated in the state."
  (list parent-state))

;; variable, state -> list of states
(defmethod generate ((variable random-variable)
		     &optional (parent-state (make-root-state variable)))
  (let* ((vertex (vertex variable))
	 (a-priori-congruent
	   (congruent-values (model variable) vertex parent-state))
	 (hidden? (not (nth-value 1 (gethash vertex parent-state))))
	 (parent-probability (gethash :probability parent-state))
	 (congruent-values
	   (if hidden? a-priori-congruent
	       (when (member (gethash vertex parent-state)
			     a-priori-congruent :test #'equal)
		 (list (gethash vertex parent-state))))))
    (flet ((normalize (probabilities)
	     (let ((total-mass (apply #'pr:add probabilities)))
	       (mapcar (lambda (p) (pr:div p total-mass)) probabilities)))
	   (branch-state (value &optional probability)
	     (format t "Parent p ~a p ~a~%" parent-probability probability) 
	     (branch-state (unless (null probability)
			     (pr:mul probability parent-probability))
			   vertex value parent-state)))
      (if *estimate?*
	  (mapcar #'branch-state congruent-values)
	  (let* ((probabilities (conditional-probabilities variable congruent-values parent-state))
		 (probabilities (if hidden? (normalize probabilities)
				    probabilities)))
	    (mapcar #'branch-state congruent-values probabilities))))))

;; model, [state] -> list of states
(defmethod generate ((m bayesian-network) &optional (parent-state (make-root-state m)))
  "Generate the congruent states that can be reached from PARENT-STATE."
  (let ((congruent-states (graph-generate m parent-state)))
    (%write-states m congruent-states)
    congruent-states))

(defmethod transition-states ((m dynamic-bayesian-network) states &optional marginal)
  ;; todo optional marginal
  (let ((marginal (union (state-variables m) marginal))) ;; or maybe this is fine
    (marginalize
     (rotate-states m (marginalize states marginal)
		    :persist marginal)
     (mapcar #'previous marginal))))

;; THere are two ways of transitioning
;; 1. Caling transition on a moment
;; 2. Calling generate on a moment, then (setf (hidden-state m) (transition-states congruent-states)

;; model, state -> list of states (side effect: update internal state)
(defmethod transition ((m dynamic-bayesian-network) observation final?
		       &optional marginal) ;; Fix optional marginal
  (let* ((*generate-a-priori-states* nil)
	 (marginal (union (state-variables m) marginal))
	 (observed-states (observe-states m (hidden-state m) observation))
	 (branched-states 
	   (mapcar (lambda (s)
		     (let ((new-states (graph-generate m s)))
		       (when final? (next-sequence m new-states))
		       (transition-states m new-states marginal)))
		   observed-states))
	 (hidden-state (apply #'append branched-states))
	 (hidden-state (marginalize hidden-state (mapcar #'previous marginal))))
    (setf (hidden-state m) hidden-state)))

;; TOD: Fix dealing with marginal

(defmethod sequence-generate ((m dynamic-bayesian-network) moment-sequence
			      &key marginal (write-header? t) (%moment 0))
  "Step through a sequence and transition the model. Throw an error if a moment
causes a dead end and no hidden states are generated.
When done, (hidden-states m) should hold the remaining congruent states
at the end of the sequence."
  (when write-header? (%write-header m))
  (when (= %moment 0) (initialize m))
  (unless (null moment-sequence)
    (let* ((*moment* %moment))
      (transition m (observation m (car moment-sequence))
		  (null (cdr moment-sequence)) marginal)
      (when (null (hidden-state m))
	(error "No states congruent with observation ~a of moment ~a at position #~a in sequence ~a." (format-obs (observation m (car moment-sequence))) (car moment-sequence)
	       *moment* *sequence*))
      (sequence-generate m (cdr moment-sequence)
			 :%moment (1+ %moment)
			 :write-header? nil
			 :marginal marginal))))

(defmethod input-generate ((m bayesian-network) input)
  (generate m (observation m input)))

(defmethod moment-generate ((m dynamic-bayesian-network) moment)
  "Generate slice of joint distribution congruent with moment."
  (flet ((branch (state)
	   (let* ((congruent-states (graph-generate m state)))
	     (if *generate-a-priori-states*
		 (congruent-states congruent-states) 
		 congruent-states))))
    (let* ((observation (observation m moment))
	   (observed-states (observe-states m (hidden-state m) observation)))
      (apply #'append (mapcar #'branch observed-states)))))

(defmethod input-generate ((m dynamic-bayesian-network) input)
  (sequence-generate m input))

(defun make-marginal-state (variables values &optional trace)
  (let ((state (make-hash-table)))
    (unless (null trace)
      (setf (gethash :trace state) trace))
    (loop for variable in variables for value in values do
      (setf (gethash variable state) value))
    state))

(defun update-marginal-state (marginal-state p &optional trace)
  (let ((state-p (state-probability marginal-state)))
    (unless (null trace)
      (let ((current-trace-p (state-probability (gethash :trace marginal-state)))
	    (new-trace-p (state-probability trace)))
      (when (> new-trace-p current-trace-p)
	(setf (gethash :trace marginal-state) trace))))
    (set-state-probability marginal-state
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
	     (probability (state-probability state))
	     (new-state (update-marginal-state marginal-state probability trace)))
	(setf (gethash values marginal) new-state)))
    (hash-table-values marginal)))

(defmethod observe-states ((m dynamic-bayesian-network) states observation)
  (mapcar (lambda (s) (observe-state m s observation)) states))

(defmethod observe-state ((m dynamic-bayesian-network) state observation)
  (maphash (lambda (vertex value) (when (member vertex (slot-value m 'vertices))
				    (setf (gethash vertex state) value)))
	   observation)
  state)

(defmethod rotate-states ((m dynamic-bayesian-network) states
			  &key
			    (keep-trace? (keep-trace? m))
			    (persist (model-variables m)))
  (unless (null states)
    (cons (rotate-state m (car states) :keep-trace? keep-trace? :persist persist)
	  (rotate-states m (cdr states) :keep-trace? keep-trace? :persist persist))))

(defmethod rotate-state ((m dynamic-bayesian-network) state
			 &key (persist (model-variables m))
			   (keep-trace? (keep-trace? m)))
  "\"Rotate\" a state. In  rotated (a priori) version of a state, every parameter
X is renamed ^X and variables of the form ^X in STATE are dropped.
Optionally, a trace can be kept which stores a hash table containing the dropped
values, their probability and the previous trace."
  (let ((new-state (make-hash-table)))
    (setf (gethash :probability new-state) (gethash :probability state))
    (when keep-trace?
      (let ((trace (make-hash-table)))
	(dolist (key (cons :probability (cons :trace persist)))
	  (setf (gethash key trace) (gethash key state)))
	(setf (gethash :trace new-state) trace)))
    (dolist (variable persist new-state)
      (setf (gethash (previous variable) new-state)
	    (gethash-or variable state)))))

(defmethod variable-value ((var random-variable) state)
  (cons (gethash (vertex var) state)
	(loop for v in (distribution-parents var) collect (gethash v state))))

(defmethod variable-values ((var random-variable) states)
  (unless (null states)
    (cons (variable-value var (car states))
	  (variable-values var (cdr states)))))

(defmethod next-sequence ((m dynamic-bayesian-network) congruent-states)
  (dolist (v (vertices m) nil)
    (next-sequence (model-variable m v) congruent-states)))

(defmethod next-sequence ((var random-variable) congruent-states)
  (let ((variable-values (variable-values var congruent-states)))
    (next-sequence (distribution var) variable-values)))

(defmethod persistent-vertices ((m dynamic-bayesian-network)
				intermediate-marginalization)
  "Persistent variables are those variables whose values are preserved between
moment. That is, they are not marginalized out."
  (let ((additional (if (null intermediate-marginalization)
			(vertices m)
			(when (listp intermediate-marginalization)
			  intermediate-marginalization))))
    (union (model-variables m) additional)))

(defun set-state-probability (state probability)
  (setf (gethash :probability state) probability))
			      
(defun state-probability (state)
  (gethash :probability state))

(defmethod evidence ((m bayesian-network) congruent-states)
  (state-probability
	   (car (marginalize congruent-states (observed-variables m)))))

(defun posterior (congruent-states)
  "Same as dividing normalizing the congruent states."
  (let* ((probabilities (mapcar #'state-probability congruent-states))
	 (evidence (apply #'+ probabilities)))
    (dolist (state congruent-states congruent-states)
      (set-state-probability state (pr:div (state-probability state) evidence)))))
			      
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

(defmethod %write-states ((m bayesian-network) states &optional (output (output m)))
  (dolist (state states)
    (%write-state m state output)))

(defmethod %write-state ((m bayesian-network) state &optional (output (output m)))
  (format output "~a,~a~{,~a~^~},~a,~a~%" *sequence* *moment*
	  (loop for v in (if (null (output-vars m))
			     (slot-value m 'vertices)
			     (output-vars m))
		collect (funcall (formatter-function m v) (gethash v state)))
	  (congruent? state) (gethash :probability state)))

(defmethod trace-back ((m dynamic-bayesian-network) state &optional (vertices (persistent-vertices m)))
  ;;  (let ((new-trace (cons (gethash variable state) trace))
  (let ((trace (gethash :trace state)))
    (if (null trace)
	state
	(let ((vertices (mapcar #'previous vertices)))
	  (dolist (v vertices)
	    (setf (gethash v state) (gethash (basename v) trace)))
	  (trace-back m state vertices)))))

(defmethod pprint-state ((m bayesian-network) state)
  (format t "(~{~{~A~^:~}~^ ~}): ~a"
	  (loop for k being each hash-key of state
		if (not (eq k :probability))
		  collect (list k (gethash k state)))
	  (gethash :probability state)))

(defmethod probability-table ((var random-variable))
  (let ((header (append (parents var)
			(list (vertex var) 'probability))))
    (cons header (probability-table (distribution var)))))

(defun state-probability-table (states &key variables sort)
  (let* ((states-variables
	   (unless (null states)
	     (loop for key being each hash-key of (car states)
		   if (not (member key '(:probability :trace)))
		     collect key)))
	 (variables (if (null variables) states-variables
			(if (subsetp variables states-variables)
			    variables
			    (progn
			      (warn "Variables not found in state: ~a"
				    (set-difference states-variables variables))
			      variables))))
	 (header (append variables (list :probability)))
	 (states (if  (not (null sort))
		      (sort states #'< :key (lambda (s) (gethash :probability s)))
		      states))
	 (table))
    (dolist (state states)
      (push (loop for key in (append variables (list :probability))
		  collect (gethash-or key state))
	    table))
    (cons header table)))
