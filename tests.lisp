(cl:defpackage #:jackdaw-tests
  (:use #:common-lisp #:jackdaw)
  (:export
   "TEST-ALL" "TEST")
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

(defmodel test-1 () ()
  ((a () (uniform ()) ())
   (b (a ^b) (uniform ()) ())
   (c () (uniform ()) ())))

(defparameter *all-tests* nil)

(defun test-function-symbol (function)
  (intern (format nil "TEST-~A"
		  (if (stringp function)
		      (string-upcase function)
		      (symbol-name function)))
	  :jackdaw-tests))

(defun test-all ()
  (dolist (test *all-tests*)
    (format t "Testing ~a... " test)
    (funcall (test-function-symbol test))
    (format t "✔~%")))

(defun test (function)
  (let ((test-function (test-function-symbol function)))
    (if (member function *all-tests*)
	(funcall test-function)
	(error "No test found for ~a" function))))

(defmacro deftest (name &body test)
  (pushnew name *all-tests*)
  `(defun ,(test-function-symbol name) ()
     ,@test))

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
    
(defmacro deftest-with-model (name (model &rest args) &body test)
  `(deftest ,name
     (let ((model (make-instance ',model ,@args)))
       ,@test)))

(deftest-with-model state-variables (test-1 :observe '(c))
  (assert (equal (jackdaw::state-variables model) '(b))))

(deftest-with-model model-variables (test-1 :observe '(c))
  (assert (equal (jackdaw::model-variables model) '(b c))))

(deftest-with-model rotate-state (test-1 :observe '(c))
  (let ((state (plist->hash-table
		'(:probability 1
		  ^a 0
		  ^b 1
		  ^c 0
		  a 2
		  b 3
		  c 4))))
    (assert (plists-equal (hash-table->plist (jackdaw::rotate-state model state :keep-trace? nil))
			  '(:probability 1
			    ^b 3
			    ^c 4)))))
					

;;(deftest-with-model estimate (test-1)
;;  (assert (equal (jackdaw::estimate model ) '(b c))))

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
    (assert (equal (length (jackdaw::marginalize states nil)) 1))
    (assert (approx-equal (pr:out (gethash :probability (car (jackdaw::marginalize states nil)))) 1))
    (assert (equal (length (jackdaw::marginalize states '(a))) 2))
    (assert (equal (gethash 'a (first (jackdaw::marginalize states '(a))))
		   0) () "Unexpected value of first marginal state. Maybe order changed?")
    (assert (approx-equal (pr:out (gethash :probability (first (jackdaw::marginalize states '(a)))))
			  .4))))

		   
			  
	   
		 
