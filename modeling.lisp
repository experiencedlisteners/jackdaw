(cl:in-package #:jackdaw)

(defvar +ngram-filler+ '✳)
(defvar +singleton+ (list +inactive+))

;; Constraint definition utility macros

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

