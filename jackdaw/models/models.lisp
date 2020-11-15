(cl:in-package #:jackdaw)

(defmodel one-var (generative-model) ()
  ((x () '(t nil))) ())

(defmodel two-var (generative-model) ()
  ((x () '(t nil))
   (y () '(t nil)))
  ())

(defmodel three-var (generative-model) ()
  ((:x () '(t nil))
   (:y () '(t nil))
   (:z (:x :y) (deterministic (and $x $y))))
  ())
   

(defmodel tonality (generative-model) ()
  ((mode (^mode) (persistent '(major minor))))
  ())

(defmodel ppm-generative-model (generative-model)
  (&key update-exclusion order-bound training? (mixtures t) (escape :c)))

(defmodel enculturation (ppm-generative-model)
  (training? ioi-domain meter-distribution)
  ((M (^m)
      (categorical () training? meter-distribution)
      (persistent (meter-domain model))
      :key (lambda (m) (list (getf m 'period) (getf m 'pulses))))
   (D (^d ^p m)
      (accumulator-model (m) training?)
      (recursive (loop for ioi in (ioi-domain model)
		       collect (cons (+ $^p ioi) $^d))
		 (deterministic '(*))))
   (P (^p m d)
      (recursive (loop for phase below (car $m)
		       collect phase)
		 (list (mod (car $d) (car $m)))))
   (IOI (d ^p)
	(recursive (list (- (car $d) $^p))
		   (list '*)))))

(defmethod meter-domain ((m enculturation))
  ;; loop over (meter-distribution m) and collect the keys
  )

(defmodel enculturation (ppm-generative-model)
  ((ioi-domain :initarg :ioi-domain :reader ioi-domain)
   (meter-distribution :initarg :meter-distribution :reader meter-distribution))
  ((M (^m)
      (categorical () (training? meter-distribution))
      (persistent (meter-domain model))
      :inputs (period pulses))
   (D (^d ^p m)
      (accumulator-model (m) (training?))
      (recursive (loop for ioi in (ioi-domain model)
		       collect (cons (+ $^p ioi) $^d))
		 (deterministic '(*))))
   (P (^p m d)
      (recursive (loop for phase below (car $m)
		       collect phase)
		 (list (mod (car $d) (car $m)))))
   (IOI (d ^p)
	(recursive (list (- (car $d) $^p))
		   (list '*))
	:inputs (ioi)))
  :required-parameters (ioi-domain meter-distribution))


(defmethod initialize-instance :after ((m enculturation)
				       update-exclusion
				       (mixtures t) order-bound (escape :c))
  (assert (not (member 0 (ioi-domain m))) nil
	  "0 may not be a member of the IOI domain.")
  (assert (or (null meter-params)
	      (< (abs (- (apply #'+ (mapcar #'cadr meter-params)) 1)) 1e-10))
	  nil "Meter probabilities do not sum to approximately one.")
  (let ((ppm-dist (distribution m 'd)))
    (setf (training? ppm-dist) (training? m)
	  (slot-value ppm-dist 'mixtures) mixtures
	  (slot-value ppm-dist 'escape) escape
	  (slot-value ppm-dist 'update-exclusion) update-exclusion
	  (slot-value ppm-dist 'order-bound) order-bound))
  (if (training? m) (funcall #'hide m) (funcall #'hide m 'M 'P))
  (warn "Training is ~A." (if (training? m) "ON" "OFF"))
  (setf (slot-value m 'meter-domain) (mapcar #'car meter-params))
  (loop for (meter probability) in meter-params do
       (set-param (distribution m 'M) nil meter probability)))


(defmodel enculturation (generative-model)
  ((ioi-domain :initarg :ioi-domain :reader ioi-domain)
   (meter-domain :initarg :meter-domain :reader meter-domain)
   (training? :initarg :training? :reader training? :initform nil))
  ((M (^m) (persistent (meter-domain model))
      :inputs (period pulses))
   (D (^d ^p)
      (recursive (loop for ioi in (ioi-domain model)
		       collect (cons (+ $^p ioi) $^d))
		 (deterministic '(*))))
   (P (^p m d)
      (recursive (loop for phase below (car $m)
		       collect phase)
		 (list (mod (car $d) (car $m)))))
   (IOI (d ^p)
	(recursive (list (- (car $d) $^p))
		   (list '*))
	:inputs (ioi)))
  ((D (m) (accumulator-model))
   (M () (categorical)))
  :required-fields (ioi-domain))

(defwriter enculturation (m)
  (list (read-from-string (with-output-to-string (s) (call-next-method m s)))
	(ioi-domain m) (meter-domain m)))
(defreader enculturation (m data)
  (let ((model (first data))
	(ioi-domain (second data))
	(meter-domain (third data)))
    (with-input-from-string (s (write-to-string model))
      (call-next-method m s)
      (setf (slot-value m 'ioi-domain) ioi-domain)
      (setf (slot-value m 'meter-domain) meter-domain))))

(defmethod initialize-instance :after ((m enculturation)
				       &key meter-params update-exclusion
				       (mixtures t) order-bound (escape :c))
  (assert (not (member 0 (ioi-domain m))) nil
	  "0 may not be a member of the IOI domain.")
  (assert (or (null meter-params)
	      (< (abs (- (apply #'+ (mapcar #'cadr meter-params)) 1)) 1e-10))
	  nil "Meter probabilities do not sum to approximately one.")
  (let ((ppm-dist (distribution m 'd)))
    (setf (training? ppm-dist) (training? m)
	  (slot-value ppm-dist 'mixtures) mixtures
	  (slot-value ppm-dist 'escape) escape
	  (slot-value ppm-dist 'update-exclusion) update-exclusion
	  (slot-value ppm-dist 'order-bound) order-bound))
  (if (training? m) (funcall #'hide m) (funcall #'hide m 'M 'P))
  (warn "Training is ~A." (if (training? m) "ON" "OFF"))
  (setf (slot-value m 'meter-domain) (mapcar #'car meter-params))
  (loop for (meter probability) in meter-params do
       (set-param (distribution m 'M) nil meter probability)))
