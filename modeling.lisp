(cl:in-package #:jackdaw)

(defvar +ngram-filler+ '✳)
(defvar +singleton+ (list +inactive+))

(defun ensure-list (variable)
  (if (inactive? variable) nil variable))

;; Constraint definition utility macros

(defmacro deterministic (congruent-value) `(list ,congruent-value))

(defmacro normal (constraint) constraint)

(defmacro recursive (^self constraint initialization-constraint)
  `(if (inactive? ,^self) ,initialization-constraint ,constraint))

(defmacro persist (^self constraint)
  `(recursive ,^self (list ,^self) ,constraint))

(defmacro one-shot (^self constraint)
  `(persist ,^self ,constraint))
 
(defmacro accumulate (^self constraint &optional initialization-constraint)
  `(recursive ,^self
    (mapcar (lambda (s) (cons s ,^self)) ,constraint)
    (mapcar #'list ,(or initialization-constraint constraint))))

(defmacro markov (order ^self constraint &optional initialization-constraint)
  `(recursive
    ,^self
    (let ((order ,(cond
		    ((numberp order)
		     `(min (length ,^self) ,order))
		    ((null order)
		     `(length ,^self))
		    (t
		     `(if (null ,order) (length ,^self)
			  (min (length ,^self) ,order))))))
      (mapcar (lambda (s) (cons s (subseq ,^self 0 order))) ,constraint))
    (mapcar #'list ,(or initialization-constraint constraint))))

(defmacro ngram (^self n constraint &optional initialization-constraint)
  `(recursive ,^self
    (mapcar (lambda (s) (cons s (subseq ,^self 0 (1- ,n)))) ,constraint)
    (mapcar (lambda (s)	(cons s (loop repeat (1- ,n) collect +ngram-filler+)))
	    ,(or initialization-constraint constraint))))

(defmacro chain (constraint &rest dependencies)
  `(if (not (every (lambda (s) (inactive? s))
		   (list ,@dependencies)))
       ,constraint
       +singleton+))
