(cl:in-package #:jackdaw)
;; Special symbols

(defvar +inactive+ '*)
(defvar +singleton+ (list +inactive+))
(defvar +ngram-filler+ '*)

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
  ((output :accessor output :initarg :output :initform nil)
   (outputvars :accessor outputvars :initform nil :initarg :outputvars)
   (variables :initarg :variables :reader variables :type hash-table
	      :initform (required-arg variables generative-model))))

(defclass random-variable ()
  ((name :initarg :name :reader name)
   (output :initarg :output :reader output)
   (hidden :initarg :hidden :accessor hidden :initform nil)
   (distribution :reader distribution)
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

;; Representation of states and variable identifiers

(defun constr-arg (v)
    "Given a keyword representation of a variable, return
its congruency constraint function symbol. For example,
given :X, return $X"
  (intern (format nil "$~A" (symbol-name v))))

(defun congruency-function (v)
  (intern (format nil "$~A" (symbol-name v)) :jackdaw))

(defun input-argument (i)
  "Given an input identifier, return
its congruency constraint function symbol. For example,
given :X, return <X"
  (intern (format nil "<~A" (symbol-name i))))

(defun getarg (key state)
  "Given a variable identifier, obtain
its value from a model state."
  (cond
    ((hash-table-p state)
     (multiple-value-bind (v found?)
	 (gethash key state)
       ;;(unless found? (warn "~A not found in state or moment." key))
       v))
    ((listp state)
     (getf state (%kw key)))))
;;(getf state key))

(defun previous (v)
  (intern (format nil "^~A" (symbol-name v))))

(defun apriori? (v)
  (eq (elt (symbol-name v) 0) #\^))

(defun basename (s)
  "Given an a priori version of a variable name, return its
stem. For example, if S is :^X, (BASENAME S) is :X."
  (intern (subseq (symbol-name s) 1) :jackdaw))

(defun get-horizontal-arguments (arglist)
  (loop for arg in arglist
     if (eq (elt (symbol-name arg) 0) #\^)
       collect (basename arg)))

(defun get-vertical-arguments (arglist)
  (loop for arg in arglist
     if (not (eq (elt (symbol-name arg) 0) #\^))
     collect arg))

;; Constraint definition macros

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

(defmacro ngram-accumulator (constraint &optional initialization-constraint n)
  `(recursive 
    (mapcar (lambda (s) (cons s (subseq $^self 0 (1- ,n)))) ,constraint)
    (mapcar (lambda (s)	(cons s (loop repeat (1- ,n) collect +ngram-filler+)))
	    ,(or initialization-constraint constraint))))

(defmacro chain (constraint dependencies)
  `(if (not (every (lambda (s) (inactive? s))
		   (list ,@(mapcar #'constr-arg dependencies))))
       ,constraint
       +singleton+))

(defmacro required-arg (symbol cls)
  `(error ,(format nil "~a is a required initialization argument of ~a."
		   (%kw symbol) cls)))
  
;; Model definition macro

(defmacro defmodel (class superclasses parameters variables)
  (let ((direct-slots (%lambda-list->direct-slots parameters class))
	(edges (make-hash-table))
	(distributions (make-hash-table))
	(var-specs) (vertices))
    (assert (or (null superclasses)
		(eq (length (remove-if (lambda (c) (subtypep c 'generative-model)) superclasses))
		    (1- (length superclasses))))
	    () "At least one superclass of ~a must be of type GENERATIVE-MODEL." class)
    `(muffle-redefinition-warnings
       (defclass ,class ,(if (null superclasses) '(generative-model) superclasses)
	 (,@direct-slots))
       ,@(loop
	   for variable in variables
	   collect
	   (destructuring-bind (v parents distribution constraint
				&key (key `(lambda (m) (getf m ',(%kw v))))
				  (formatter '#'identity) hidden)
	       variable
	     (push v vertices)
	     (setf (gethash v distributions) distribution)
	     (setf (gethash v edges) parents)
	     (push (list v key formatter hidden) var-specs)
	     `(progn
		(defmethod ,(congruency-function v) ((model ,class) args)
		  (declare (ignorable model))
		  (multiple-value-bind (,@(mapcar #'constr-arg parents))
		      (apply #'values args)
		    (declare (ignorable ,@(mapcar #'constr-arg parents)))
		    (handler-case
			,(if (member (previous v) parents)
			     `(let (($^self ,(constr-arg (previous v))))
				(declare (ignorable $^self))
				,constraint)
			     constraint)
		      (error (e)
			(warn "Error in a priori congruency constraint of ~A!" ',v)
			(error e))))))))
       (defmethod vertices ((m ,class))
	 ',vertices)
       (defmethod edge-table ((m ,class)) ,edges)
       (defun ,(intern (format nil "MAKE-~A-INSTANCE" (symbol-name class)))
	   (,@parameters ,@(unless (member '&key parameters) '(&key)) output outputvars)
	 ,(let ((x (gensym)))
	    `(let ((,x (make-hash-table)))
	       ,@(loop
		   for (v key formatter hidden) in var-specs
		   collect
		   `(setf (gethash ',v ,x)
			  (make-instance
			   'jackdaw::random-variable
			   :hidden ,hidden
			   :name ',v :key ,key :output ,formatter)))
	       (let ((,class (make-instance ',class :variables ,x
					    :output output :outputvars outputvars
					    ,@(%lambda-list->plist parameters))))
		 ,@(loop
		     for v in vertices
		     collect
		     (destructuring-bind (dist dist-args &rest dist-params)
			 (gethash v distributions)
		       (assert (subsetp dist-args (gethash v edges))
			       () "Arguments to distribution of ~a must be subset of parents." v)
		       (assert (subtypep dist 'distribution) ()
			       "Distribution argument of the variable ~a (~a) must be of type DISTRIBUTION." 
			       v dist)
		       `(let ((distribution
				(make-instance ',dist :arguments ',dist-args
						      :variable ',v
						      ,@dist-params)))
			  (setf (slot-value (gethash ',v (variables ,class)) 'distribution)
				distribution))))
		 ,class)))))))


;; Model mechanics

(defmethod formatter-function ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((v (gethash variable (variables m))))
    (output v)))

(defmethod variable-constraint ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((v (gethash variable (variables m))))
    (constraint v)))

(defmethod observed-value ((m generative-model) variable moment)
  "Obtain the observed value of VARIABLE from MOMENT given M."
  (let ((v (gethash variable (variables m))))
    (funcall (key v) moment)))

(defmethod hidden? ((m generative-model) variable)
  (hidden (gethash variable (variables m))))

(defmethod observed? ((m generative-model) v)
  (not (hidden? m v)))

(defmethod hide ((m generative-model) &rest variables)
  "Hide VARIABLES."
  (dolist (v variables)
    (setf (hidden (gethash v (variables m))) t)))

(defmethod make-observable ((m generative-model) &rest variables)
  "Make VARIABLES observable."
  (dolist (v variables)
    (setf (hidden (gethash v (variables m))) nil)))

(defmethod get-var-distribution ((m generative-model) variable)
  (distribution (gethash variable (variables m))))

(defmethod marginal-params ((m generative-model))
  (let* ((horizontal-dependencies
	 (mapcar (lambda (v) (get-horizontal-arguments (edges m v))) (vertices m))))
    (remove-duplicates (apply #'append horizontal-dependencies))))

(defmethod a-posteriori-congruent ((m generative-model) variable parents-state moment state)
  "Return the a posteriori congruent states of VARIABLE."
  (if (hidden? m variable) t
      (equal state (observed-value m variable moment))))
;;  (let ((constraint (posterior-constraint m variable)))
;;    (format t "Var: ~A. Constraint: ~A~%" variable constraint)
;;    (if (not (or (null constraint) (hidden? m variable)))
;;	(handler-case
;;	    (funcall constraint m parents-state moment state)
;;	  (error (e)
;;	    (warn "Error in a priori congruency constraint of ~A!" variable)
;;	    (error e)))
;;	t)))

(defmethod congruent-states ((m generative-model) variable parents-state)
  (let ((parents (edges m variable)))
    (funcall (congruency-function variable) m
	     (loop for p in parents collect (gethash p parents-state)))))

(defmethod generate-states ((m generative-model) variables previous-state moment
			    &optional (final? t))
  (let* ((root? (null (cdr variables)))
	 (parent-states (if root? (list previous-state)
			    (generate-states m (cdr variables) previous-state moment nil)))
	 (variable (car variables))
	 (new-states))
    ;;(warn "Generating ~a.~%" variable)
    (dolist (parents-state parent-states new-states)
      (let* ((probability (gethash :probability parents-state))
	     (congruent-states (congruent-states m variable parents-state))
	     (distribution (probabilities (get-var-distribution m variable)
					  parents-state congruent-states))
	     (parent-congruent? (or root? (gethash :congruent? parents-state))))
	;;(warn "root? ~s parent congruent? ~s p=~a" root? parent-congruent? probability)
	(when (eq (length congruent-states) 0)
	  (warn "~A has no a priori congruent states." variable))
	(dolist (s congruent-states)
	  (let ((congruent? (and parent-congruent?
				 (a-posteriori-congruent m variable parents-state moment s)))
		(new-state (copy-hash-table parents-state))
		(s-probability (pr:mul probability (gethash s distribution))))
	    (setf (gethash variable new-state) s)
	    (setf (gethash :congruent? new-state) congruent?)
	    ;;(warn "~a ~a ~a ~a" final? s probability s-probability)
	    (when final? (write-state m new-state congruent? s-probability))
	    ;; Incongruent final states should not be returned,
	    ;; otherwise, do generate them.
	    (when (or (not final?) congruent?)
	      (when (null s-probability)
		(warn "State ~a not found in distribution." s))
	      (setf (gethash :probability new-state) s-probability)
	      (push new-state new-states))))))))

(defmethod rotate-state ((m generative-model) state &key (keep-trace? t))
  "\"Rotate\" a state. In  rotated (a priori) version of a state, every parameter
:X is renamed :^X and variables of the form :^X in STATE are dropped."
  (let ((new-state (make-hash-table)))
    (setf (gethash :probability new-state) (gethash :probability state))
    (when keep-trace?
      (let ((trace (make-hash-table)))
	(dolist (key (cons :trace (marginal-params m)))
	  (setf (gethash key trace) (gethash key state)))
	(setf (gethash :trace new-state) trace)))
    (dolist (variable (marginal-params m) new-state)
      (setf (gethash (previous variable) new-state) (gethash variable state)))))

(defmethod write-header ((m generative-model) &optional (output (output m)))
  (format output "sequence,moment~{,~a~^~},congruent,probability~%" 
	  (loop for v in (if (null (outputvars m))
			     (vertices m)
			     (outputvars m)) collect (string-downcase (symbol-name v)))))

(defmethod write-state ((m generative-model) state congruent probability
			&optional (output (output m)))
  (format output "~a,~a~{,~a~^~},~a,~a~%" *sequence* *moment*
	  (loop for v in (if (null (outputvars m))
			     (vertices m)
			     (outputvars m)) collect (funcall (formatter-function m v)
							      (getarg v state)))
	  congruent probability))

(defun trace-back (state variable &optional trace)
  (let ((new-trace (cons (gethash variable state) trace))
	(previous-state (gethash :trace state)))
    (if (null previous-state)
	new-trace
	(trace-back previous-state variable new-trace))))

(defmethod root-state ((m generative-model))
  "The root state. Every root nodes in a Bayesian networks is conditioned
on this root state."
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) (pr:in 1))
    (dolist (variable (mapcar #'previous (marginal-params m)) state)
      (setf (gethash variable state) +inactive+))))

(defun marginalize (states variables &optional (marginal (make-hash-table :test #'equal)))
  "Given a list of states, create a hash table representing a marginal distribution."
  (dolist (state states marginal)
    (let* ((state-prob (gethash :probability state))
	   (trace (gethash :trace state))
	   (key (loop for v in variables collect (gethash v state)))
	   (marginal-prob (gethash key marginal))
	   (new-prob (if (null marginal-prob) state-prob
			 (pr:add state-prob (car marginal-prob)))))
      (setf (gethash key marginal) (cons new-prob trace)))))

(defun marginal->states (marginal variables)
  (let ((states))
    (maphash (lambda (state-key p)
	       (let ((state (make-hash-table)))
		 (setf (gethash :probability state) (car p))
		 (setf (gethash :trace state) (cdr p))
		 (loop for v in variables
		    for s in state-key do
		      (setf (gethash v state) s))
		 (push state states)))
	     marginal)
    states))

(defmethod process-dataset ((m generative-model) dataset &optional (write-header? t))
  (let* ((sequence (car dataset))
	 (*sequence* (car sequence)))
    ;;(warn "~a" (car sequence))
    (process-sequence m (cdr sequence) write-header?))
  (unless (null (cdr dataset))
    (process-dataset m (cdr dataset) nil)))
	 
(defmethod process-sequence ((m generative-model) moments
			     &optional
			       (write-header? t)
			       (moment 0)
			       (congruent-states (list (root-state m))))
  (when write-header? (write-header m))
  (if (null moments)
      (dolist (v (vertices m) congruent-states)
	(next-sequence (get-var-distribution m v) congruent-states))
      (let* ((*moment* moment)
	     (new-congruent-states (transition m (car moments) congruent-states)))
	(process-sequence m (cdr moments) nil (1+ moment) new-congruent-states))))

(defmethod transition ((m generative-model) moment congruent-states)
  (let ((marginal (make-hash-table :test #'equal))
	(marginal-variables (mapcar #'previous (marginal-params m))))
    (dolist (previous-state congruent-states)
      (let* ((new-states (model-congruent-states m previous-state moment))
	     (marginal-states (mapcar (lambda (s) (rotate-state m s)) new-states)))
	(marginalize marginal-states marginal-variables marginal)))
    (marginal->states marginal marginal-variables)))

(defmethod model-congruent-states ((m generative-model) previous-state moment)
  "Call GENERATE-STATES on topologically sorted list of variables, from which
inactive variables have been pruned."
  (let ((a-posteriori-congruent-states
	  (generate-states m (get-vertical-arguments
			      (topological-sort m)) previous-state
			      moment)))
    (when (eq (length a-posteriori-congruent-states) 0)
      (warn "Set of a posteriori congruent states is empty at moment
~a at position #~a in sequence ~a." moment *moment* *sequence*))
    a-posteriori-congruent-states))

(defmethod list-congruent-states ((m generative-model) previous-state moment)
  (let ((states (model-congruent-states m previous-state moment)))
    (loop for state in states
	  collect
	  (append (loop for v in (vertices m) collect (gethash v state))
		  (list (gethash :probability state))))))

(defmethod pprint-state ((m generative-model) state)
  (format t "(~{~{~A~^:~}~^ ~}): ~a~%" 
	  (loop for v in (vertices m)
		collect (list v (gethash (previous v) state)))
	  (gethash :probability state)))

(defmethod evidence ((m generative-model) congruent-states)
  (apply #'probabilities:add
	 (loop for s in congruent-states collect (gethash :probability s))))
  
;; Model serialization

(defwriter generative-model (m)
    (let ((params))
      (loop for var being the hash-key using (hash-value dist) of (variables m) do
	(setf (getf params var) (serialize dist)))
      params))
(defreader generative-model (m distributions)
  (loop for (var dist) on distributions by #'cddr do
       (with-input-from-string (s (write-to-string dist))
	 (deserialize (getf (distributions m) var) s))))

;; Optimizations:
;; Sort vertices after model creation
;; Represents states as vectors of length |VERTICES|
;; Create ONE hash table that returns for each vertex the index of
;; its state in a state-vector
