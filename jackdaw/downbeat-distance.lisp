(cl:in-package #:jackdaw)

(defmodel downbeat-distance (generative-model)
  ((ioi-domain :initarg :ioi-domain :reader ioi-domain)
   (meter-domain :reader meter-domain)
   (training? :initarg :training? :reader training? :initform nil))
  ((M (^m)
      (one-shot (meter-domain model))
      :inputs (period pulses) :posterior-constraint
      (recursive t (equal $m (list <period <pulses))))
   (B (^b ^p m)
      (chain
       (accumulator
	(loop for ioi in (ioi-domain model) collect (+ $^p ioi)))
       (^p))
      :output (lambda (s) (if (eq +inactive+ s) "" (car s)))
      :inputs (ioi)
      :posterior-constraint
      (chain-posterior (eq (+ $^p <ioi) (car $b)) (^p)))
   (period (m) (deterministic (car $m)))
   (pulses (m) (deterministic (cadr $m)))
   (P (^p m b) ;; Observe this to constrain the first phase.
      (recursive (deterministic (mod (car $b) (car $m)))
		 (loop for phase below (car $m) collect phase))
      :inputs (ioi)
      :posterior-constraint
      (recursive t (eq $p <ioi)))
   (P0 (^p0 p) (one-shot (deterministic $p)))
   (IOI (b ^p) (chain (deterministic (- (car $b) $^p)) (^p))))
  ((B (m) (accumulator-model))
   (M nil (categorical)))
  :required-fields (ioi-domain))


(defwriter downbeat-distance (m)
  (list (read-from-string (with-output-to-string (s) (call-next-method m s)))
	(ioi-domain m) (meter-domain m)))
(defreader downbeat-distance (m data)
  (let ((model (first data))
	(ioi-domain (second data))
	(meter-domain (third data)))
    (with-input-from-string (s (write-to-string model))
      (call-next-method m s)
      (setf (slot-value m 'ioi-domain) ioi-domain)
      (setf (slot-value m 'meter-domain) meter-domain))))

(defmethod initialize-instance :after ((m downbeat-distance)
				       &key meter-params update-exclusion
					 (mixtures t) order-bound (escape :c))
  (assert (not (member 0 (ioi-domain m))) nil
	  "0 may not be a member of the IOI domain.")
  (assert (or (null meter-params)
	      (< (abs (- (apply #'+ (mapcar #'cadr meter-params)) 1)) 1e-10))
	  nil "Meter probabilities do not sum to approximately one.")
  (let ((ppm-dist (distribution m 'b)))
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

;; test reading and writing
;; why is generate state ran more than once? or is it?

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
