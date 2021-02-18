(cl:in-package #:jackdaw)
;; Special symbols

(defvar +inactive+ '✳)
(defvar +singleton+ (list +inactive+))
(defvar +ngram-filler+ '✳)

;; Globals for enabling meta-programming with models

(defparameter *models* nil) ; list of defined models
(defparameter *model-parameters* nil) ; plist of parameters per model

;; Settings

(defparameter *generate-a-priori-states* t
  "When T, only all a-priori congruent states are generated.
This is useful when we want to calculate things like the entropy of
the predictive distribution.")

(defun model-exists? (model-symbol-or-name &optional (package :jackdaw))
  (not (null (find-model model-symbol-or-name package))))

(defun find-model (model-symbol-or-name &optional (package :jackdaw))
  (let ((symbol (find-symbol (if (stringp model-symbol-or-name)
				 (string-upcase model-symbol-or-name)
				 (symbol-name model-symbol-or-name))
			     package)))
    (when (subtypep symbol 'generative-model)
      symbol)))

;; Global parameters for CSV output writing.

(defparameter *sequence* nil)
(defparameter *moment* nil)

;; Object serialization helper macros.

(defmacro defreader (type (model data) &body body)
  `(defmethod deserialize ((,model ,type) stream)
     (let ((,data (read stream)))
       ,@body)))

(defmacro defwriter (type (model) data)
  `(defmethod serialize ((,model ,type) &optional stream)
     (let ((data ,data))
       (unless (null stream) (write data :stream stream))
       data)))

;; Objects

(defclass generative-model (dag)
  ((%var-specs :allocation :class :type list)
   (%parameter-slots :reader %parameter-slots :allocation :class :type list)
   ;;(%dist-specs :allocation :class :type list)
   (output :accessor output :initarg :output :initform nil)
   (output-vars :accessor output-vars :initform nil :initarg :output-vars)
   (distributions :reader distributions :type hash-table)
   (variables :reader variables :type hash-table)))

(defclass random-variable ()
  ((name :initarg :name :reader name)
   (output :initarg :output :reader output)
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

(defun congruency-function (v)
  (intern (format nil "$~A" (symbol-name v)) :jackdaw))

(defun previous (v)
  (intern (format nil "^~A" (symbol-name v)) (symbol-package v)))

(defun horizontal? (v)
  (eq (elt (symbol-name v) 0) #\^))

(defun basename (s)
  "Given an a priori version of a variable name, return its
stem. For example, if S is :^X, (BASENAME S) is :X."
  (intern (subseq (symbol-name s) 1) (symbol-package s)))

(defmethod horizontal-edges ((m generative-model) vertex)
  (loop for v in (edges m vertex)
     if (horizontal? v)
       collect (basename v)))

(defmethod vertical-edges ((m generative-model) vertex)
  (loop for v in (edges m vertex)
     if (horizontal? v)
       collect (basename v)))

(defun get-vertical-arguments (vertices)
  (loop for v in vertices
     if (not (horizontal? v))
     collect v))

;; Constraint definition utility macros

(defun inactive? (s)
  (eq s +inactive+))

(defmacro deterministic (congruent-value) `(list ,congruent-value))

(defmacro normal (constraint) constraint)

(defmacro recursive (constraint initialization-constraint)
  `(if (inactive? $^self) ,initialization-constraint ,constraint))

(defmacro persistent (constraint)
  `(recursive (list $^self) ,constraint))

(defmacro one-shot (constraint)
  `(persistent ,constraint))

(defmacro accumulator (constraint &optional initialization-constraint)
  `(recursive
    (mapcar (lambda (s) (cons s $^self)) ,constraint)
    (mapcar #'list ,(or initialization-constraint constraint))))

(defmacro ngram (constraint n &optional initialization-constraint)
  `(recursive 
    (mapcar (lambda (s) (cons s (subseq $^self 0 (1- ,n)))) ,constraint)
    (mapcar (lambda (s)	(cons s (loop repeat (1- ,n) collect +ngram-filler+)))
	    ,(or initialization-constraint constraint))))

(defmacro chain (constraint dependencies)
  `(if (not (every (lambda (s) (inactive? s))
		   (list ,@(mapcar #'constr-arg dependencies))))
       ,constraint
       +singleton+))

;; Model definition macro

(defmacro required-arg (symbol cls)
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
  (let* ((direct-slots (%lambda-list->direct-slots parameters class))
	 (parameter-names (mapcar #'%param-name (remove '&key parameters)))
	 (edges (make-hash-table))
	 (dist-specs) (var-specs) (vertices)
	 (methods))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c) (subtypep c 'generative-model)) superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type GENERATIVE-MODEL." class)
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
	 `(defmethod ,(congruency-function v) ((model ,class) args)
	    (declare (ignorable model))
	    (destructuring-bind (,@(mapcar #'constr-arg parents))
		args
	      (declare (ignorable ,@(mapcar #'constr-arg parents)))
	      (let (,@(loop for parameter in parameter-names collect
			    `(,parameter (slot-value model ',parameter)))
		    ,@(when (member (previous v) parents)
			`(($^self ,(constr-arg (previous v))))))
		(declare (ignorable
			  ,@parameter-names
			  ,@(when (member (previous v) parents) `($^self))))
		(handler-case
		    ,constraint
		  (error (e)
		    (warn "Error in a priori congruency constraint of ~A!" ',v)
		    (error e))))))
	 methods)))
    `(muffle-redefinition-warnings
       (pushnew (%kw ',class) *models*)
       (setf (getf *model-parameters* ',class) ',parameters)
       (defclass ,class ,(if (null superclasses) '(generative-model) superclasses)
	 ((%var-specs :initform ',var-specs)
	  (%parameter-slots :initform ',parameter-names)
	  ;;(%dist-specs :initform ',dist-specs)
	  ,@direct-slots))
       ,@methods
       (defmethod vertices ((m ,class))
	 ',vertices)
       (defmethod edge-table ((m ,class)) ,edges)
       (defmethod initialize-instance :after ((model ,class) &key)
	 (let ((distributions (make-hash-table))
		,@(loop for p in (mapcar #'%param-name (remove '&key parameters))
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
		      `(setf (gethash ',v distributions)
			     (apply #'make-instance ',dist
				    ,(append `(list :arguments ',dist-args
						    :variable-symbol ',v)
					     dist-params)))))
	    (setf (slot-value model 'distributions) distributions)))
       (defun ,(intern (format nil "MAKE-~A-MODEL" (symbol-name class)))
	   (,@parameters ,@(unless (member '&key parameters) '(&key)) output output-vars observe)
	 (make-instance ',class :output output :output-vars output-vars :observe observe
			,@(%lambda-list->plist parameters))))))

(defmethod initialize-instance :after ((model generative-model) &key observe)
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
		    :hidden hidden
		    :name v :key (eval key) :output (eval formatter))))
    (setf (slot-value model 'variables) variables))
  (apply #'observe (cons model observe)))

;; Model serialization

(defwriter generative-model (m)
	   (let ((dist-params)
		 (parameter-slots (loop for s in (%parameter-slots m)
					collect (slot-value m s))))
	     (loop for var being the hash-key using (hash-value dist) of (distributions m) do
	       (setf (getf dist-params var) (serialize dist)))
	     (list parameter-slots dist-params)))
(defreader generative-model (m data)
  (destructuring-bind (parameters dist-params)
      data
    (loop for p in parameters
	  for s in (%parameter-slots m)
	  do (setf (slot-value m s) p))
    (loop for (var dist) on dist-params by #'cddr do
      (with-input-from-string (s (write-to-string dist))
	(deserialize (gethash var (distributions m)) s)))))

;; Jackdaw mechanics

(defmethod formatter-function ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((v (gethash variable (variables m))))
    (output v)))

(defmethod observed-value ((m generative-model) variable moment)
  "Obtain the observed value of VARIABLE from MOMENT given M."
  (let ((v (gethash variable (variables m))))
    ;; TODO wrap in handler case and warn about errors
    (funcall (key v) moment)))

(defmethod hidden? ((m generative-model) variable)
  (hidden (gethash variable (variables m))))

(defmethod %set-hidden ((m generative-model) hidden vertices)
  (dolist (v vertices)
    (let ((variable
	    (gethash-or v (variables m)
		       (error "~a is not a variable of ~a" v (type-of m)))))
      (setf (hidden variable) hidden))))

(defmethod hide ((m generative-model) &rest vertices)
  "Hide VARIABLES."
  (%set-hidden m t vertices))
  
(defmethod observe ((m generative-model) &rest vertices)
  "Make VARIABLES observable."
  (%set-hidden m nil vertices))

(defmethod get-var-distribution ((m generative-model) variable)
  (gethash variable (distributions m)))

(defmethod state-variables ((m generative-model))
  (reduce #'union (mapcar (lambda (v) (horizontal-edges m v)) (vertices m))))

(defmethod model-variables ((m generative-model))
  (union (state-variables m) (observed-variables m)))

(defun get-probability (state distribution)
  (let ((probability (gethash-or state distribution
			 (warn "State ~a not found in distribution." state))))
    probability))

(defmethod order ((m generative-model))
  (length (vertices m)))

(defmethod observations ((m generative-model) moment)
  (let ((observations (make-hash-table)))
    (dolist (v (observed-variables m) observations)
      (setf (gethash v observations) (observed-value m v moment)))))

(defmethod observed-states (states observations)
  "Return the subset of STATES that is consistent with OBSERVATIONS."
  (loop for s in states if (observed? s observations) collect s))

(defmethod observed? (state observations)
  "Return T if STATE is consistent with OBSERVATIONS."
  (every #'identity
	 (loop for v being each hash-key of observations
	       collect (equal (gethash v observations) (gethash-or v state)))))

(defmethod observed-variables ((m generative-model))
  (remove-if (lambda (v) (hidden? m v)) (vertices m)))

(defmethod congruent-variable-states ((m generative-model) variable parents-state)
  "Apply a variable's congruency function to its dependencies to obtain the
congruent states of the variable."
  (let* ((parents (edges m variable))	 
	 (arguments
	   (loop for p in parents
		 collect (gethash-or p parents-state)))
	 (states (funcall (congruency-function variable) m arguments)))
    (when (eq (length states) 0)
      (warn "~A has no a priori congruent states." variable))
    states))

(defmethod make-root-state ((m generative-model))
  "Every root node in a Bayesian network is implicitly conditioned
on this root state."
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) (pr:in 1))
    (dolist (variable (state-variables m))
      (setf (gethash variable state) +inactive+))
    state))

(defun make-state (probability &optional old-state)
  (let ((new-state (if (null old-state) (make-hash-table)
		       (copy-hash-table old-state))))
    (setf (gethash :probability new-state) probability)
    new-state))

(defparameter *estimate?* nil)

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
		 (states (generate-moment m observation (rotate-state m state))))
	    (assert (> (length states) 0) ()
		    "Could not generate observation ~a of moment ~w" (format-obs observation) moment)
	    (assert (eq (length states) 1) ()
		    "Cannot estimate from model that is not fully observed.")
	    	    (push (loop for v in (vertices m) collect (gethash v (car states)))
		  observation-sequence)
	    (setf state (car states))))
	(push (reverse observation-sequence) observation-dataset)))
    (reverse observation-dataset)))

(defmethod estimate ((m generative-model) dataset)
  (warn "generating observations")
  (let* ((dataset (generate-observations m dataset)))
    (dotimes (vertex-index (order m) m)
      (let* ((v (elt (vertices m) vertex-index))
	     (distribution (get-var-distribution m v))
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

(defmethod probability ((m generative-model) observation)
  (evidence m (generate-sequence m observation)))
	      
"(defmethod probabilities ((m generative-model) states)
  (let ((distributions (make-hash-table :test #'equal))
	(probabilities))
    (flet ((probability (var val state)
	     (let ((distribution (get-var-distribution m var))
		   (arguments (mapcar (lambda (arg) (gethash arg state))
				      (arguments distribution))))
	       (multiple-value-bind (distribution found?)
		   (gethash (cons var arguments) distributions)
		 (unless found?
		   (setf (gethash (cons var arguments) distributions)
			 (probabilities distribution arguments 
    (dolist (state states)
      (dolist (v (vertices m))
	(
"
		     
(defmethod generate-moment ((m generative-model) observations
			    &optional (state (make-root-state m))
			      (variables (get-vertical-arguments (topological-sort m))))
  "It is possible to separate the generation of states from the generation of probabilities.
However, that would require generating the distributions while generating states (as below)
and keeping track of them in a hash table that is returned by generate moment. These
distributions can then be used to look up probabilities in a separate function."
  (let* ((parent-states (if (null (cdr variables)) (list state)
			    (generate-moment m observations state (cdr variables))))
	 (variable (car variables))
	 (new-states))
    ;;(format t "Generating ~a~%" variable)
    (multiple-value-bind (value observed?)
	(gethash variable observations)
      (dolist (parents-state parent-states new-states)
	(let* ((a-priori-congruent-states
		 (congruent-variable-states m variable parents-state))
	       (constrain-states?
		 (and observed? (not *generate-a-priori-states*)))
	       (congruent-states
		 (if constrain-states?
		     (when (member value a-priori-congruent-states)
		       (list value))
		     a-priori-congruent-states))
	       (probability (unless *estimate?* (gethash :probability parents-state)))
	       (distribution (unless *estimate?*
			       (probabilities (get-var-distribution m variable)
					      parents-state a-priori-congruent-states))))
	  (dolist (s congruent-states)
	    ;;(when (eq variable 'letter-observer)
	    ;;  (format t "~a prob: ~a X ~a = ~a~%" variable probability (get-probability s distribution)
	;;	      (pr:mul probability (get-probability s distribution))))
	    (let ((new-state (make-state (unless *estimate?*
					   (pr:mul probability (get-probability s distribution)))
					 parents-state)))
	      (setf (gethash variable new-state) s)
	      (push new-state new-states))))))))

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

(defmethod rotate-state ((m generative-model) state &key (keep-trace? t))
  "\"Rotate\" a state. In  rotated (a priori) version of a state, every parameter
X is renamed ^X and variables of the form ^X in STATE are dropped."
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
  
(defmethod transition ((m generative-model) moment congruent-states &optional keep)
  (let ((observations (observations m moment))
	(persistent-variables (union (model-variables m) keep)))
    (flet ((generate-moment (state)
	     (let* ((states (generate-moment m observations
					     (rotate-state m state)))
		    (congruent-states
		      (if *generate-a-priori-states*
			  (observed-states states observations)
			  states)))
	       (when (eq (length congruent-states) 0)
		 (warn "No a posteriori congruent states at moment
~a at position #~a in sequence ~a." moment *moment* *sequence*))
	       (write-states m congruent-states observations)
	       (marginalize congruent-states persistent-variables))))
      (let* ((new-states
	       (apply #'append (mapcar #'generate-moment congruent-states))))
	(marginalize new-states persistent-variables)))))

(defmethod generate-dataset ((m generative-model) dataset &optional (write-header? t))
  (let* ((sequence (car dataset))
	 (*sequence* (car sequence)))
    ;;(warn "~a" (car sequence))
    (generate-sequence m (cdr sequence) write-header?))
  (unless (null (cdr dataset))
    (generate-dataset m (cdr dataset) nil)))
	 	 
(defmethod generate-sequence ((m generative-model) moments
			      &key
				keep
			       (write-header? t)
			       (moment 0)
			       (congruent-states (list (make-root-state m))))
  (when write-header? (write-header m))
  (if (null moments)
      (dolist (v (vertices m) congruent-states)
	(next-sequence (get-var-distribution m v) congruent-states))
      (let* ((*moment* moment)
	     (new-congruent-states (transition m (car moments) congruent-states keep)))
	;;(format t "Evidence ~a, prob first con state: ~a n-cong: ~a~%"
	;;	(evidence m new-congruent-states)
	;;	(gethash :probability (car new-congruent-states))
	;;	(length new-congruent-states))
	(generate-sequence m (cdr moments)
			   :keep keep
			   :write-header? nil
			   :moment (1+ moment)
			   :congruent-states new-congruent-states))))

(defmethod evidence ((m generative-model) congruent-states)
  (gethash :probability
	   (car (marginalize congruent-states nil))))

;;(defmethod posterior-distribution ((m generative-model) congruent-states variables)
;;  (let ((states (marginalize congruent-states (mapcar #'previous variables)))
;;	(evidence (evidence m congruent-states)))
;;    (loop for param being the hash-key using (hash-value prob) of table do
;;      (setf (gethash param table) (/ (car prob) evidence)))
;;    table))

(defmethod write-header ((m generative-model) &optional (output (output m)))
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

(defmethod write-states ((m generative-model) states observations &optional (output (output m)))
  (dolist (state states)
    (write-state m state observations output)))

(defmethod write-state ((m generative-model) state observations &optional (output (output m)))
  (format output "~a,~a~{,~a~^~},~a,~a~%" *sequence* *moment*
	  (loop for v in (if (null (output-vars m))
			     (vertices m)
			     (output-vars m))
		collect (funcall (formatter-function m v) (gethash v state)))
	  (observed? state observations) (gethash :probability state)))

(defun trace-back (state variable &optional trace)
  (let ((new-trace (cons (gethash variable state) trace))
	(previous-state (gethash :trace state)))
    (if (null previous-state)
	new-trace
	(trace-back previous-state variable new-trace))))

(defmethod format-hash-table (hash-table)
  (format nil "(~{~{~A~^:~}~^ ~})"
	  (loop for k being each hash-key of hash-table
		collect (list k (gethash k hash-table)))))

(defmethod pprint-state ((m generative-model) state)
  (format t "(~{~{~A~^:~}~^ ~})"
	  (loop for k being each hash-key of state
		collect (list k (gethash k state))))
  (format t "(~{~{~A~^:~}~^ ~}): ~a~%" 
	  (append (loop for v in (vertices m)
			collect (list v (gethash v state)))
		  (loop for v in (vertices m)
			collect (list (previous v) (gethash (previous v) state))))
	  (gethash :probability state))
  state)

;; Optimizations:
;; Sort vertices after model creation
;; Represents states as vectors of length |VERTICES| and represent variables as indices
;; Create ONE hash table that returns for each vertex the index of
;; its state in a state-vector
