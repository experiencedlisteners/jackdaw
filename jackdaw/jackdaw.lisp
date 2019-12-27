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

(defvar +inactive+ '*)
(defvar +singleton+ (list +inactive+))
(defvar +ngram-filler+ '*)

;; Global parameters for CSV output writing.

(defparameter *sequence* nil)
(defparameter *event* nil)

;; Utility

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

(defun hash-table->alist (hashtab)
  (let ((alist '()))
    (maphash #'(lambda (key value) (setq alist (cons (list key value) alist)))
             hashtab)
    alist))

(defun alist->hash-table (alist &key (test #'equal))
  (let ((hashtable (make-hash-table :test test)))
    (mapc #'(lambda (x) (setf (gethash (car x) hashtable) (cadr x)))
          alist)
    hashtable))

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
   (hidden :accessor hidden :initform nil)
   (required-fields :reader required-fields :initform nil)
   (observed :reader observed :initform nil :type 'list)
   (output-functions :reader output-functions :initform nil :type 'list)
   (prior-constraints :reader prior-constraints :type 'hashtable)
   (posterior-constraints :reader posterior-constraints :type 'hashtable)
   (distributions :accessor distributions :type 'list)))

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

(defun constraint-argument (v)
    "Given a keyword representation of a variable, return
its congruency constraint function symbol. For example,
given :X, return $X"
    (intern (format nil "$~A" (symbol-name v))))

(defun input-argument (i)
  "Given an input identifier, return
its congruency constraint function symbol. For example,
given :X, return <X"
  (intern (format nil "<~A" (symbol-name i))))

(defun tokeyword (symbol)
  (intern (symbol-name symbol) :keyword))

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
     (getf state (tokeyword key)))))
;;(getf state key))

(defun apriori (v)
  (intern (format nil "^~A" (symbol-name v)) :jackdaw))

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

(defmacro make-constraint (self parents constraint &optional inputs)
  `(lambda (model %state &optional moment ,(constraint-argument self))
     (declare (ignorable model moment ,(constraint-argument self)))
     (let* (,@(loop for v in parents collect
		   (list (constraint-argument v) `(getarg ',v %state)))
	    ,@(when (member (apriori self) parents)
		    `(($^self ,(constraint-argument (apriori self)))))
	      ,@(loop for i in inputs collect
		     (list (input-argument i) `(getarg ',i moment))))
	      ,constraint)))

(defmacro deterministic (congruent-value) `(list ,congruent-value))

(defmacro normal (constraint) constraint)

(defmacro recursive (constraint initialization-constraint)
  `(if (eq $^self +inactive+) ,initialization-constraint ,constraint))

(defmacro one-shot (constraint)
  `(recursive (list $^self) ,constraint))

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
  `(if (not (every (lambda (s) (eq s +inactive+))
		   (list ,@(mapcar #'constraint-argument dependencies))))
       ,constraint
       +singleton+))

(defmacro chain-posterior (constraint dependencies)
  `(if (not (every (lambda (s) (eq s +inactive+))
		   (list ,@(mapcar #'constraint-argument dependencies))))
       ,constraint
       t))


;; Model definition macro

(defmacro defmodel (class superclasses direct-slots variable-specs distribution-specs
		    &key required-fields)
  (let* ((edges (make-hash-table))
	 (variables) (posterior-constraints) (prior-constraints)
	 (output-functions) (distributions))
    (loop for specification in variable-specs do
	 (destructuring-bind
	       (v parents prior-constraint &key posterior-constraint inputs
		  (output '#'identity))
	     specification
	   (push v variables)
	   (push output output-functions)
	   (setf (gethash v edges) parents)
	   (let ((prior-constraint (macroexpand `,prior-constraint))
		 (posterior-constraint (macroexpand `,posterior-constraint)))
	     (push (macroexpand `(make-constraint ,v ,parents ,prior-constraint))
		   prior-constraints)
	     (push
	      (unless (null posterior-constraint)
		(macroexpand `(make-constraint ,v ,parents ,posterior-constraint ,inputs)))
	      posterior-constraints))))
    (loop for specification in distribution-specs do
	 (destructuring-bind (v parents (type &rest args))
	     specification
	   (setf (gethash v edges)
		 (union (gethash v edges) parents))
	   (push `(',v (make-instance ',type :arguments '(,@parents)
				      :variable ',v ,@args))
		 distributions)))
    `(defclass ,class ,superclasses
       ((distributions :initform (list ,@(apply #'append (reverse distributions))))
	(vertices :initform '(,@(reverse variables)))
	(edge-table :initform ,edges)
	(output-functions :initform (list ,@(reverse output-functions)))
	(prior-constraints :initform (list ,@(reverse prior-constraints)))
	(posterior-constraints :initform (list ,@(reverse posterior-constraints)))
	(required-fields :initform '(,@required-fields))
	,@direct-slots))))

;; Model mechanics

(defmethod initialize-instance :after ((m generative-model) &key)
  (when (not (every (lambda (v) (member v (vertices m))) (outputvars m)))
    (error "Some outputvars are not a variable in the model."))
  (when (not (every (lambda (s) (slot-boundp m s)) (required-fields m)))
    (error "~A requires ~{~a~^,~} to be provided as initialization arguments."
	   (type-of m) (required-fields m))))

(defmethod output-function ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((p (position variable (vertices m))))
    (elt (output-functions m) p)))

(defmethod prior-constraint ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((p (position variable (vertices m))))
    (elt (prior-constraints m) p)))

(defmethod posterior-constraint ((m generative-model) variable)
  "Return the congruency constraint associated with VARIABLE in model M."
  (let ((p (position variable (vertices m))))
    (elt (posterior-constraints m) p)))

(defmethod hidden? ((m generative-model) variable)
  (member variable (hidden m)))

(defmethod hide ((m generative-model) &rest variables)
  "Hide VARIABLES. Make VARIABLES the only hidden variables. 
To observe everything, call without variables."
  (loop for v in variables if (not (member v (vertices m))) do
       (warn "Hiding a variable, ~a, that is not part of the ~a model."
	     v (type-of m)))
  (setf (hidden m) variables))

(defmethod distribution ((m generative-model) variable)
  (getf (distributions m) variable (make-instance 'uniform)))

(defmethod marginal-params ((m generative-model))
  (let* ((horizontal-dependencies
	 (mapcar (lambda (v) (get-horizontal-arguments (edges m v))) (vertices m))))
    (remove-duplicates (apply #'append horizontal-dependencies))))

(defmethod a-priori-congruent ((m generative-model) variable parents-state)
    (handler-case
	(funcall (prior-constraint m variable) m parents-state)
      (error (e)
	(warn "Error in a priori congruency constraint of ~A!" variable)
	(error e))))

(defmethod observed? ((m generative-model) v)
  (not (or (null (posterior-constraint m v)) (hidden? m v))))

(defmethod a-posteriori-congruent ((m generative-model) variable parents-state moment state)
  "Return the a posteriori congruent states of VARIABLE."
  (let ((constraint (posterior-constraint m variable)))
    (if (not (or (null constraint) (hidden? m variable)))
	(handler-case
	    (funcall constraint m parents-state moment state)
	  (error (e)
	    (warn "Error in a priori congruency constraint of ~A!" variable)
	    (error e)))
	t)))

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
	     (a-priori (a-priori-congruent m variable parents-state))
	     (distribution (probabilities (distribution m variable)
					  parents-state a-priori))
	     (parent-congruent? (or root? (gethash :congruent? parents-state))))
	;;(warn "root? ~s parent congruent? ~s p=~a" root? parent-congruent? probability)
	(when (eq (length a-priori) 0)
	  (warn "~A has no a priori congruent states." variable))
	(dolist (s a-priori)
	  (let ((congruent? (and parent-congruent?
				 (a-posteriori-congruent m variable parents-state moment s)))
		(new-state (copy-hash-table parents-state))
		(s-probability (pr:mul probability (gethash s distribution))))
	    (setf (gethash variable new-state) s)
	    (setf (gethash :congruent? new-state) congruent?)
	    ;;(warn "~a ~a ~a" s probability s-probability)
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
      (setf (gethash (apriori variable) new-state) (gethash variable state)))))

(defmethod write-header ((m generative-model) &optional (output (output m)))
  (format output "sequence,event~{,~a~^~},congruent,probability~%" 
	  (loop for v in (if (null (outputvars m))
			     (vertices m)
			     (outputvars m)) collect (string-downcase (symbol-name v)))))

(defmethod write-state ((m generative-model) state congruent probability
			&optional (output (output m)))
  (format output "~a,~a~{,~a~^~},~a,~a~%" *sequence* *event*
	  (loop for v in (if (null (outputvars m))
			     (vertices m)
			     (outputvars m)) collect (funcall (output-function m v)
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
    (dolist (variable (mapcar #'apriori (marginal-params m)) state)
      (setf (gethash variable state) +inactive+))))

(defun marginalize (states variables &optional (marginal (make-hash-table :test #'equal)))
  "Given a list of states, create a hash table representing a marginal distribution."
  (dolist (state states)
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
			       (event 0)
			       (congruent-states (list (root-state m))))
  (when write-header? (write-header m))
  (if (null moments)
      (dolist (v (vertices m) congruent-states)
	(next-sequence (distribution m v) congruent-states))
      (let* ((*event* event)
	     (new-congruent-states (transition m (car moments) congruent-states)))
	(process-sequence m (cdr moments) nil (1+ event) new-congruent-states))))

(defmethod transition ((m generative-model) moment congruent-states)
  (let ((marginal (make-hash-table :test #'equal))
	(marginal-variables (mapcar #'apriori (marginal-params m))))
    (dolist (previous-state congruent-states)
      (let* ((new-states (model-congruent-states m previous-state moment))
	     (marginal-states (mapcar (lambda (s) (rotate-state m s)) new-states)))
	(marginalize marginal-states marginal-variables marginal)))
    (marginal->states marginal marginal-variables)))

(defmethod model-congruent-states ((m generative-model) previous-state moment)
  "Call GENERATE-STATES on topologically sorted list of variables, from which
inactive variables have been pruned."
  (let ((congruent-states
	 (generate-states m (get-vertical-arguments
			     (topological-sort m)) previous-state
			     moment)))
    (when (eq (length congruent-states) 0)
      (warn "Posterior congruency constraints not satisfied in any state for input
~a at moment #~a." moment *event*))
    congruent-states))

;; Model serialization

(defwriter generative-model (m)
    (let ((params))
      (loop for (var dist) on (distributions m) by #'cddr collect
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
