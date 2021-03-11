(cl:defpackage #:jackdaw-tests
  (:use #:common-lisp #:jackdaw)
  (:export
   "RUN-ALL" "RUN")
  (:documentation "Unittests for jackdaw."))

(in-package #:jackdaw-tests)

;;(jackdaw:defmodel test-model () ()
;;  ((A () (jackdaw:bernouilli () :p (/ 3 4)) (list nil t))))

;; With no observations this should just generate two states with probabilities 1/4 and 3/4
;; With observations this should generate whatever state is observed
;; With *generate-a-priori-states* the probabilities should be the same
;; but a priori states should be printed

;;(jackdaw:defmodel test-model-2 () ()
;;  ((A (^a) (jackdaw:bernouilli () :p (/ 3 4)))))

;; With observations (nil t) this should generate a state with symbol t and probability 3/16
;; 

(defparameter *all-tests* nil)

(defun test-function-symbol (function)
  (intern (format nil "TEST-~A"
		  (if (stringp function)
		      (string-upcase function)
		      (symbol-name function)))
	  :jackdaw-tests))

(defun add-test (name test)
  (setf (getf *all-tests* (intern (symbol-name name) :keyword))
	test))

(defun find-test (symbol-or-name)
  (let ((symbol (if (stringp symbol-or-name)
		    (intern (string-upcase symbol-or-name) :keyword)
		    (intern (symbol-name symbol-or-name) :keyword))))
    (if (null (getf *all-tests* symbol))
	(error "No test called ~a appears to exist" symbol-or-name)
	(getf *all-tests* symbol))))

(defun run-all ()
  (loop for (name) on *all-tests* by #'cddr collect 
    (run name)))

(defun run (test)
  (format t "Testing ~a " test)
  (funcall (find-test test))
  (format t " ✔~%"))

(defmacro deftest (name &body test)
  `(add-test ',name (lambda () (progn ,@test))))
    
(defmacro deftest-with-model (name (model &rest args) &body test)
  `(deftest ,name
     (let ((,model (make-instance ',model ,@args)))
       ,@test)))

(defmacro deftest-with-models (name (&rest models) &body test)
  `(deftest ,name
     (let (,@(loop for (model args) in models
		   collect `(,model (make-instance ',model ,@args))))
       ,@test)))

(defmacro test (test &optional message arguments)
  `(progn
     (assert ,test () ,message ,arguments)
     (format t ".")))

(defun plist->hash-table (plist &key (test #'eql))
  (let ((hash-table (make-hash-table :test test)))
    (loop for (key value) on plist by #'cddr do
      (setf (gethash key hash-table) value))
    hash-table))

(defun hash-table->plist (hash-table)
  (let ((plist))
    (loop for key being each hash-key of hash-table do
      (push (gethash key hash-table) plist)
      (push key plist))
    plist))

(defun plist-keys (plist)
  (loop for key in plist by #'cdr collect key))

(defun sets-equal (set-a set-b &key (test #'eql))
  (and (eq (length set-a) (length set-b))
       (eq (length (intersection set-a set-b :test test))
	   (length set-a))))

(defun plists-equal (plist-a plist-b &key (key-test #'eql) (value-test #'eql))
  (and (sets-equal (plist-keys plist-a) (plist-keys plist-b)
		   :test key-test)
       (every (lambda (k) (funcall value-test
				   (getf plist-a k)
				   (getf plist-b k)))
	      (plist-keys plist-a))))

(defun hash-tables-equal (hash-table-a hash-table-b &key (key-test #'eql) (value-test #'eql))
  (plists-equal (hash-table->plist hash-table-a)
		(hash-table->plist hash-table-b)
		:key-test key-test
		:value-test value-test))

(defun states-equal (a b)
  (let ((a (copy-list a))
	(b (copy-list b))
	(p-a (getf a :probability))
	(p-b (getf b :probability)))
    (remf a :probability)
    (remf b :probability)
    (and (approx-equal p-a p-b)
	 (plists-equal a b))))

(defun state-lists-equal (a b)
  (sets-equal (mapcar #'hash-table->plist a)
	      (mapcar #'hash-table->plist b)
	      :test #'states-equal))

(defmacro test-state-lists-equal (a b)
  `(let ((a (mapcar #'hash-table->plist ,a))
	 (b (mapcar #'hash-table->plist ,b)))
     (test (sets-equal a b :test #'states-equal))))

(defmacro test-hash-tables-equal (a b &key (key-test '#'eql)
					(value-test '#'eql))
  `(let ((a (hash-table->plist ,a))
	 (b (hash-table->plist ,b)))
     (test (plists-equal a b :key-test ,key-test :value-test ,value-test))))


;;;;;;;; Unit tests ;;;;;;;;;;;

;; Bayesian networks

;; (defnetwork bn-1 () () (:a (:b :a)))

;; (set-distribution bn-1 :a ()
;; 		  (make-bernouilli-distribution :p :psymbol))
;; (set-distribution bn-1 :b (:a)
;; 		  (make-cpt-distribution
;; 		   :alist-cpt '(((:p :x) . .8)
;; 				((:q :x) . .2)
;; 				((:p :y) . .3)
;; 				((:q :y) . .7))))

;; (set-values bn-1 :a

;; (set-observer bn-1 :a #'first)
;; (set-observer bn-1 :b #'second)

;; TODO: Test both with *generate-a-priori* and without

(defmodel bn-1 () ()
  ((:a () (bernouilli () :p 0.8 :psymbol :x)
       '(:x :y)
    :observer #'first)
   (:b (:a) (cpt (:a) 
		 :alist-cpt '(((:p :x) . .8)
		   ((:q :x) . .2)
		   ((:p :y) . .3)
		   ((:q :y) . .7)))       
       '(:p :q)
    :observer #'second)))

(deftest-with-model bayesnet-joint-probability (bn-1)
  (test (approx-equal (probability bn-1 ()) 1))
  (observe bn-1)
  (test (approx-equal (probability bn-1 '(:x :p))
		      (* .8 .8)))
  (test (approx-equal (probability bn-1 '(:x :q))
		      (* .8 .2)))
  (test (approx-equal (probability bn-1 '(:y :p))
		      (* .2 .3)))
  (test (approx-equal (probability bn-1 '(:y :q))
		      (* .2 .7))))

(deftest-with-model bayesnet-marginal-probability (bn-1)
  (observe bn-1 :b)
  (test (approx-equal (probability bn-1 '(nil :p))
		      (+ (* .8 .8) (* .2 .3))))
  (test (approx-equal (probability bn-1 '(nil :q))
		      (+ (* .8 .2) (* .2 .7))))
  (hide bn-1)
  (observe bn-1 :a)
  (test (approx-equal (probability bn-1 '(:y)) 0.2))
  (test (approx-equal (probability bn-1 '(:x)) 0.8)))

(deftest-with-model bayesnet-posterior-probability (bn-1)
  (observe bn-1 :b)
  (test-state-lists-equal
   (posterior (generate bn-1 (observation bn-1 '(nil :p)))) ; P(A | B = p)
   (list (plist->hash-table
	  (list :b :p :a :x
		:probability (/ (* .8 .8) (+ (* .8 .8) (* .2 .3)))))
	 (plist->hash-table
	  (list :b :p :a :y
		:probability (/ (* .2 .3) (+ (* .8 .8) (* .2 .3)))))))
  (test-state-lists-equal
   (posterior (generate bn-1 (observation bn-1 '(nil :q)))) ; P(A | B = q)
   (list (plist->hash-table
	  (list :b :q :a :x
		:probability (/ (* .8 .2) (+ (* .8 .2) (* .2 .7)))))
	 (plist->hash-table
	  (list :b :q :a :y
		:probability (/ (* .2 .7) (+ (* .8 .2) (* .2 .7)))))))
  (hide bn-1) (observe bn-1 :a)
  (test-state-lists-equal
   (posterior (generate bn-1 (observation bn-1 '(:x nil)))) ; P(B | A = x)
   (list (plist->hash-table (list :a :x :b :p :probability .8))
	 (plist->hash-table (list :a :x :b :q :probability .2))))
    (test-state-lists-equal
   (posterior (generate bn-1 (observation bn-1 '(:y nil)))) ; P(B | A = y)
   (list (plist->hash-table (list :a :y :b :p :probability .3))
	 (plist->hash-table (list :a :y :b :q :probability .7)))))

;; Dynamic Bayesian Networks

;; DBN version of BN-1
(defmodel dbn-1 (dynamic-bayesian-network) ()
  ((:a () (bernouilli () :p 0.8 :psymbol :x)
       '(:x :y)
    :observer #'first)
   (:b (:a) (cpt (:a) 
		 :alist-cpt '(((:p :x) . .8)
		   ((:q :x) . .2)
		   ((:p :y) . .3)
		   ((:q :y) . .7)))       
       '(:p :q)
    :observer #'second)))

(deftest-with-models dbn-joint-probability ((dbn-1) (bn-1))
  ;; Tests with a fully observed model
  (test (approx-equal (probability dbn-1 '(())) 1))
  (observe dbn-1)
  (observe bn-1)
  (test (approx-equal (probability dbn-1 '((:x :p)))
   		      (probability bn-1 '(:x :p))))
  (test (approx-equal (probability dbn-1 '((:x :p) (:y :q) (:x :q)))
		      (pr:mul (probability bn-1 '(:x :p))
			      (probability bn-1 '(:y :q))
			      (probability bn-1 '(:x :q))))))

(deftest-with-models dbn-marginal-probability ((dbn-1) (bn-1))
  ;; Tests with a partially observed model
  (observe dbn-1 :b)
  (observe bn-1 :b)
  (test (approx-equal (probability dbn-1 '((nil :p)))
		      (probability bn-1 '(nil :p))))
  (test (approx-equal (probability dbn-1 '((nil :q) (nil :p)))
		      (pr:mul (probability bn-1 '(nil :q))
			      (probability bn-1 '(nil :p))))))

(defmodel dbn-2 (dynamic-bayesian-network) ()
  ((:a () (uniform ())
       '(:x :y)
    :observer #'first)
   (:b (:a :^a) (uniform (:a))
       (list (list :p $a $^a) (list :q $a $^a)))
   (:c (:b) (uniform ())
       (list (car $b))
    :observer (lambda (m) (if (listp m) (second m) m)))))

(deftest-with-model dbn-vertices (dbn-2)
  (test (sets-equal (jackdaw::vertices dbn-2) '(:a :b :c :^a))))

(deftest-with-model dbn-congruent-values (dbn-2)
  (test (sets-equal (jackdaw::congruent-values
		     dbn-2 :a (make-hash-table))
		    '(:x :y) :test #'equal))
  (test (sets-equal (jackdaw::congruent-values
		     dbn-2 :b (plist->hash-table '(:a :x :^a :y)))
		    '((:p :x :y) (:q :x :y)) :test #'equal))
  (test (sets-equal (jackdaw::congruent-values
		     dbn-2 :c (plist->hash-table '(:b (:p :x :y))))
		    '(:p) :test #'equal)))

(deftest-with-model dbn-generate (dbn-2 :observe '(:c))
  (let ((posterior-congruent-states (generate dbn-2 (observation dbn-2 :p))))
    (test (equal (length posterior-congruent-states) 2))))

(deftest-with-model state-variables (dbn-2)
  (test (equal (jackdaw::state-variables dbn-2) '(:a))))

(deftest-with-model model-variables (dbn-2)
  (test (equal (jackdaw::model-variables dbn-2) '(:a)))
  (observe dbn-2 :c)
  (test (equal (jackdaw::model-variables dbn-2) '(:a :c))))

(deftest-with-model rotate-state (dbn-2)
  (let ((state (plist->hash-table
		'(:probability 1 :^a 0 :^b 1 :^c 0 :a 2 :b 3 :c 4))))
    (test (plists-equal (hash-table->plist (jackdaw::rotate-state dbn-2 state :keep-trace? nil))
			'(:probability 1 :^a 2)))
    (observe dbn-2 :b)
    (test (plists-equal (hash-table->plist (jackdaw::rotate-state dbn-2 state :keep-trace? nil))
			'(:probability 1 :^a 2 :^b 3)))))

(defmodel dbn-3 (dynamic-bayesian-network) ()
  ((:a () (bernouilli ())
       '(:x :y)
    :observer #'first)
   (:b (:a) (cpt (:a))
       '(:p :q)
    :observer (lambda (m) (if (listp m) (second m) m)))))

(deftest-with-model estimate (dbn-3 :observe '(:a :b))
  (jackdaw::estimate dbn-3 '(((:x :p) (:x :q) (:y :p) (:y :q) (:y :q))))
  (test (approx-equal (pr:out (probability (jackdaw::model-variable-distribution dbn-3 :a) '(:y)))
		      (/ 3 5)))
  (test (approx-equal (pr:out (probability (jackdaw::model-variable-distribution dbn-3 :b) '(:q :y)))
		      (/ 2 3))))

(deftest variable-value
  (let ((variable (make-instance 'random-variable :vertex 'test :distribution-parents '(a b))))
    (test (equal (jackdaw::variable-value variable (plist->hash-table '(test 1 a 2 b 3)))
		   '(1 2 3)))))

(deftest variable-values
  (let ((variable (make-instance 'random-variable :distribution-parents '(a b) :vertex 'test)))
    (test (equal (jackdaw::variable-values
		    variable
		    (list (plist->hash-table '(test x a y b z))
			  (plist->hash-table '(test p a q b r))
			  (plist->hash-table '(test u a v b w))))
		   '((x y z)
		     (p q r)
		     (u v w))))))

(defmodel dbn-4 (dynamic-bayesian-network) ()
  ((:a (:^a) (uniform ())
       '(x y)
       :observer #'first)
   (:b () (uniform ())
       '(p q)
    :observer #'second)))

(deftest-with-model transition (dbn-4)
  (let ((states
	  (list (plist->hash-table '(:a x :b p :probability .1))
		(plist->hash-table '(:a y :b q :probability .2))))
	(states-2
	  (list (plist->hash-table '(:a x :b p :trace 'root :probability .2))
		(plist->hash-table '(:a y :b q :trace 'root :probability .1))))
	(obs (observation dbn-4 nil))
	(transitioned-states
	  (list (plist->hash-table '(:a x :b p :probability .075))
		(plist->hash-table '(:a y :b p :probability .075))
		(plist->hash-table '(:a x :b q :probability .075))
		(plist->hash-table '(:a y :b q :probability .075)))))
    ;; 
    ;; Transition from STATES without trace.
    (setf (jd::hidden-state dbn-4) states)
    (setf (jd::keep-trace? dbn-4) nil)
    (test-state-lists-equal
     (jd:transition dbn-4 obs nil)
     transitioned-states)
    ;; Transition from STATES with intermediate marginalization
    (setf (jd::hidden-state dbn-4) states)
    (setf (jd::intermediate-marginalization dbn-4) t)
    (test-state-lists-equal
     (jd:transition dbn-4 obs nil)
     (list (plist->hash-table '(:probability .3))))
    ;; Transition from STATES with specific intermediate marginalization
    (setf (jd::hidden-state dbn-4) states)
    (setf (jd::intermediate-marginalization dbn-4) '(:b))
    (test-state-lists-equal
     (jd:transition dbn-4 obs states)
     (list (plist->hash-table '(:b p :probability .15))
	   (plist->hash-table '(:b q :probability .15))))
    ;; Transition from STATES with trace and specific intermediate-marginalization
    (setf (jd::hidden-state dbn-4) states)
    (setf (jd::intermediate-marginalization dbn-4) '(:a))
    (setf (jd::keep-trace? dbn-4) t)
    (test-hash-tables-equal
     (gethash :trace (car (jd:transition dbn-4 obs nil)))
     (plist->hash-table '(:probability .2 :a y :trace nil)))
    ;; Using STATES-2, test that trace contains max prob state
        (setf (jd::hidden-state dbn-4) states-2)
    (setf (jd::intermediate-marginalization dbn-4) '(:a))
    (test-hash-tables-equal
     (gethash :trace (car (jd:transition dbn-4 obs nil)))
     (plist->hash-table '(:probability .2 :a x :trace 'root)))))

(deftest cpt-alist-parameters
  (let ((d (make-instance 'cpt 
			  :alist-cpt
			  (list (cons '(happy win)  0.8)
				(cons '(sad win)    .2)
				(cons '(happy lose) 0.4)
				(cons '(sad lose)   0.6)))))
    (assert (sets-equal (domain d) '(happy sad)))
    (assert (eq (probability d '(happy win)) 0.8))
    (assert (eq (probability d '(sad lose)) 0.6))))

(defdistribution test-distribution ()
    (&key (a 'x) (b 'y) (c 'z)) ())

(deftest serialize-distribution
  (let* ((d (make-instance 'test-distribution))
	 (data (jackdaw::serialize d)))
    (assert (equal data '(x y z)))
    (let ((d2 (make-instance 'test-distribution)))
      (with-input-from-string (s (format nil "~A" data))
	(jackdaw::deserialize d2 s))
      (with-slots (a b c) d2
	  (assert (equal (list a b c) '(x y z)))))))

;;(deftest-with-model estimate (dbn-2)
;;  (test (equal (jackdaw::estimate model ) '(b c))))

(defun make-state (probability variables values)
  (let ((state (make-hash-table)))
    (setf (gethash :probability state) probability)
    (loop for var in variables for v in values do
      (setf (gethash var state) v))
    state))

(defun approx-equal (a b &optional (tolerance 1e-7))
  (< (abs (- a b)) tolerance))

(deftest marginalize
  (let ((states
	  (list
	   (make-state (pr:in .1) '(a b) '(0 0))
	   (make-state (pr:in .3) '(a b) '(0 1))
	   (make-state (pr:in .4) '(a b) '(1 0))
	   (make-state (pr:in .2) '(a b) '(1 1)))))
    (test (equal (length (jackdaw::marginalize states nil)) 1))
    (test (approx-equal (pr:out (gethash :probability (car (jackdaw::marginalize states nil)))) 1))
    (test (equal (length (jackdaw::marginalize states '(a))) 2))
    (test (equal (gethash 'a (first (jackdaw::marginalize states '(a))))
		   0) () "Unexpected value of first marginal state. Maybe order changed?")
    (test (approx-equal (pr:out (gethash :probability (first (jackdaw::marginalize states '(a)))))
			  .4))))
(deftest ngram-model
  (let ((dist (make-ngram-model-distribution))
	(data '(((a a)) ((b a)) ((a a)) ((a b)) ((b b)) ((a b))))
	(c-data `(((a a) x) ((b a) x) ((b a) x) ((a a) y) ((b a) y))))
    (estimate dist data)
    (test (approx-equal (conditional-probability dist '(a a)) 2/3))
    (test (approx-equal (conditional-probability dist '(b a)) 1/3))
    (test (approx-equal (conditional-probability dist '(a b)) 2/3))
    (test (approx-equal (conditional-probability dist '(b b)) 1/3))
    (estimate dist c-data)
    (test (approx-equal (conditional-probability dist '(a a) '(x)) 1/3))
    (test (approx-equal (conditional-probability dist '(b a) '(x)) 2/3))
    (test (approx-equal (conditional-probability dist '(a a) '(y)) 1/2))
    (test (approx-equal (conditional-probability dist '(b a) '(y)) 1/2))))
      
      
    
    
