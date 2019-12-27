(cl:in-package #:jackdaw)

(defmodel viewpoint (generative-model)
  ((basic-domain :reader basic-domain :initarg :basic-domain)
   (basic-feature :reader basic-feature)
   (training? :initarg :training? :reader training? :initform nil))
  () () :required-fields (basic-domain))

(defconstant +undefined+ '@ "The undefined symbol.")

(defmacro viewpoint-accumulator (constraint &optional initialization-constraint)
  `(recursive
    (mapcar (lambda (s) (if (eq (car $^self) +undefined+)
			    (cons s (cdr $^self))
			    (cons s $^self)))
	    ,constraint)
    (mapcar (lambda (s) (list s))
	    ,(or initialization-constraint constraint))))

(defmethod initialize-instance :after ((m viewpoint)
				       &key (mixtures t) update-exclusion
					 (escape :c) order-bound)
  (let ((ppm-dist (distribution m 'v)))
    (setf (training? ppm-dist) (training? m)
	  (slot-value ppm-dist 'mixtures) mixtures
	  (slot-value ppm-dist 'escape) escape
	  (slot-value ppm-dist 'update-exclusion) update-exclusion
	  (slot-value ppm-dist 'order-bound) order-bound)))

(defwriter viewpoint (m)
  (list (read-from-string (with-output-to-string (s) (call-next-method m s)))
	(basic-domain m)))
(defreader viewpoint (m data)
  (with-input-from-string (s (write-to-string (first data)))
    (call-next-method m s)
    (setf (slot-value m 'basic-domain) (second data))))

(defclass viewpoint-model (accumulator-model) ())

(defmethod probabilities ((d viewpoint-model) parents-state congruent-states)
  "Obtain the location object of the appropriate PPM model given context.
Context is obtained by accessing the previous self of the current variable, 
which if the current variable is an accumulator, must represent the context.
Note that PARENTS-STATE represents a state in the current moment in which any
parent variables are instantiated."
  (if (or (equal congruent-states (list +inactive+))
	  (equal congruent-states (list nil)))
      (let ((table (make-hash-table :test #'equal)))
	(setf (gethash (car congruent-states) table) (pr:in 1))
	table)
      (call-next-method)))

(defmacro defviewpoint (name input viewpoint)
  `(defmodel ,name (viewpoint)
     ((basic-feature :initform '(,input)))
     ((e (^e v)
	 (accumulator
	  (loop for e in (basic-domain model) if
	       (equal
		(funcall ,viewpoint
			 (cons e (if (equal $^e +inactive+) nil $^e)))
		(car $v))
	     collect e))
	 :inputs (,input)
	 :posterior-constraint
	 (equal ,(input-argument input) (car $e)))
      (v (^v ^e)
	 (viewpoint-accumulator
	  (remove-duplicates
	   (loop for e in (basic-domain model) collect
		(funcall ,viewpoint (cons e (if (equal $^e +inactive+) nil $^e))))
	   :test #'equal))
	 :output #'car))
     ((v () (viewpoint-model)))))

(defviewpoint ioi-vp ioi
  (lambda (events) (car events)))

(defviewpoint ppm symbol
    (lambda (events) (car events)))

(defclass ioi-vp-ignore-first (ioi-vp) ())

(defmethod process-sequence ((m ioi-vp-ignore-first) moments
			     &optional
			       write-header?
			       (event 0)
			       (congruent-states (list (root-state m))))
  (if (eq event 0)
      (process-sequence m (cdr moments) write-header? 1 congruent-states)
      (call-next-method)))

(defviewpoint delta-ioi-vp ioi
  (lambda (events) (if (null (cdr events))
		       +undefined+
		       (- (car events) (cadr events)))))

(defviewpoint ioi-ratio ioi
  (lambda (events) (- (car events) (cadr events))))

(defviewpoint ioi-contour ioi
  (lambda (events) (- (car events) (cadr events))))

(defviewpoint onset-vp ioi
  (lambda (events) (apply #'+ events)))

(defmodel ioi (viewpoint)
  ((viewpoint :initform 'ioi))
  ((IOI (^ioi) (accumulator (ioi-domain model)) :inputs (ioi)
	:posterior-constraint (eq <ioi (car $ioi))))
  ((IOI () (accumulator-model))))

(defmodel d-ioi (ioi-based)
  ((viewpoint :initform 'd-ioi))
  ((IOI (^ioi d-ioi)
	(recursive (list (+ $^ioi (car $d-ioi)))
		   (ioi-domain model))
	:inputs (ioi)
	:posterior-constraint
	(normal (eq $ioi <ioi)))
   (D-IOI (^d-ioi ^ioi)
	  (chain (accumulator (loop for ioi in (ioi-domain model) collect
				   (- ioi $^ioi))) (^ioi))))
  ((D-IOI () (accumulator-model))))

(defmodel ioi-ratio (ioi-based)
  () ((IOI (^ioi d-ioi)
	   (recursive (list (* $^ioi (car $d-ioi)))
		      (ioi-domain model))
	   :inputs (ioi)
	   :posterior-constraint
	   (normal (eq $ioi <ioi)))
      (IOI-RATIO (^ioi-ratio ^ioi)
		 (chain (accumulator (loop for ioi in (ioi-domain model)
					collect (/ ioi $^ioi))) (^ioi))))
  ((IOI-RATIO () (accumulator-model))))


(defmodel ioi-contour (ioi-based)
  () ((IOI (^ioi contour)
	   (recursive (loop for ioi in (ioi-domain model)
			 if (funcall (case (car $contour) (-1 #'<) (0 #'=) (1 #'>))
				     ioi $^ioi) collect ioi)
		      (ioi-domain model))
	   :inputs (ioi)
	   :posterior-constraint
	   (normal (eq <ioi $ioi)))
      (CONTOUR (^contour ^ioi)
	       (chain (accumulator
		       (loop for ioi in (ioi-domain model) collect
			    (cond ((> ioi $^ioi) 1)
				  ((< ioi $^ioi) -1)
				  ((= ioi $^ioi) 0))))
		      (^ioi))))
  ((IOI-RATIO () (accumulator-model))))

(defmodel onset (ioi-based)
  () ((ONSET (^onset)
	     (accumulator
	      (recursive (loop for ioi in (ioi-domain model) collect (+ (car $^onset) ioi))
			 (ioi-domain model)))
	     :inputs (ioi)
	     :posterior-constraint
	     (recursive (eq (+ (car $^onset) <ioi) (car $onset))
			(eq <ioi (car $onset)))))
  ((ONSET () (accumulator-model))))


(defmodel linked-viewpoint (ioi-based)
  ((constituents :reader constituents :initarg :constituents))
  ()
  ())
