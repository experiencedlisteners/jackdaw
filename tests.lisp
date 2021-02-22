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
     (let ((model (make-instance ',model ,@args)))
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


;;;;;;;; Unit tests ;;;;;;;;;;;

(defmodel test-1 (dynamic-bayesian-network) ()
  ((:a () (uniform ())
       '(:x :y)
    :key #'first)
   (:b (:a :^a) (uniform (:a))
       (list (list :p $a $^a) (list :q $a $^a)))
   (:c (:b) (uniform ())
      (list (car $b))
    :key (lambda (m) (if (listp m) (second m) m)))))

(deftest-with-model congruent-values (test-1)
  (test (sets-equal (jackdaw::congruent-values
		     model :a (make-hash-table))
		    '(:x :y) :test #'equal))
  (test (sets-equal (jackdaw::congruent-values
		     model :b (plist->hash-table '(:a :x :^a :y)))
		    '((:p :x :y) (:q :x :y)) :test #'equal))
  (test (sets-equal (jackdaw::congruent-values
		     model :c (plist->hash-table '(:b (:p :x :y))))
		    '(:p) :test #'equal)))

(deftest-with-model generate (test-1 :observe '(:c))
  (let ((posterior-congruent-states (generate model '(:p))))
    (test (equal (length posterior-congruent-states) 2))))

(deftest-with-model state-variables (test-1 :observe '(:c))
  (test (equal (jackdaw::state-variables model) '(:a))))

(deftest-with-model model-variables (test-1 :observe '(:c))
  (test (equal (jackdaw::model-variables model) '(:a :c))))

(deftest-with-model rotate-state (test-1 :observe '(:c))
  (let ((state (plist->hash-table
		'(:probability 1
		  :^a 0
		  :^b 1
		  :^c 0
		  :a 2
		  :b 3
		  :c 4))))
    (test (plists-equal (hash-table->plist (jackdaw::rotate-state model state :keep-trace? nil))
			  '(:probability 1
			    :^a 2
			    :^c 4)))))

(defmodel test-2 (dynamic-bayesian-network) ()
  ((:a () (bernouilli ())
       '(:x :y)
    :key #'first)
   (:b (:a) (cpt (:a))
       '(:p :q)
    :key (lambda (m) (if (listp m) (second m) m)))))

(deftest-with-model estimate (test-2 :observe '(:a :b))
  (jackdaw::estimate model '(((:x :p) (:x :q) (:y :p) (:y :q) (:y :q))))
  (test (approx-equal (exp (probability (jackdaw::model-variable-distribution model :a) '(:y)))
		      (/ 3 5)))
  (test (approx-equal (exp (probability (jackdaw::model-variable-distribution model :b) '(:q :y)))
		      (/ 2 3))))

;;(deftest-with-model estimate (test-1)
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

		   
			  
	   
		 
