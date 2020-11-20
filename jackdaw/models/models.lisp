(cl:in-package #:jackdaw)

;(defmodel one-var (generative-model) ()
;  ((x () '(t nil))) ())

;(defmodel two-var (generative-model) ()
;  ((x () '(t nil))
;   (y (x) '(t nil)))
;  ())

;(defmodel three-var (generative-model) ()
;  ((:x () '(t nil))
;   (:y () '(t nil))
;   (:z (:x :y) (deterministic (and $x $y))))
;  ())
   

(defmodel ppm-generative-model (generative-model)
  (&key update-exclusion order-bound training? (mixtures t) (escape :c)) ())

(defmodel enculturation (ppm-generative-model)
  (training? ioi-domain meter-distribution) ;; TODO: identify parameters PER distribution
  ((M (^m)
      (categorical ()
		   :training? (training? enculturation)
		   :p (meter-distribution enculturation))
      (persistent (meter-domain model))
      :key (lambda (m) (list (getf m 'period) (getf m 'pulses))))
   (D (^d ^p m)
      (accumulator-model (m)
			 :training? (training? enculturation))
      (recursive (loop for ioi in (ioi-domain model)
		       collect (cons (+ $^p ioi) $^d))
		 (deterministic '(*))))
   (P (^p m d)
      (uniform ())
      (recursive (loop for phase below (car $m)
		       collect phase)
		 (list (mod (car $d) (car $m)))))
   (IOI (d ^p ^ioi)
	(uniform ())
	(recursive (list (- (car $d) $^p))
		   (list '*)))))

(defmethod meter-domain ((m enculturation))
  ;; loop over (meter-distribution m) and collect the keys
  )

(defmethod initialize-instance :after ((m enculturation) &key)
  (assert (not (member 0 (ioi-domain m))) nil
	  "0 may not be a member of the IOI domain.")
  (assert (or (null (meter-distribution m))
	      (< (abs (- (apply #'+ (mapcar #'cadr (meter-distribution m))) 1)) 1e-10))
	  nil "Meter probabilities do not sum to approximately one.")
  (if (training? m) (funcall #'hide m) (funcall #'hide m 'M 'P))
  (warn "Training is ~A." (if (training? m) "ON" "OFF")))

;; TODO: Integrate this!
;;  (loop for (meter probability) in meter-params do
;;       (set-param (distribution m 'M) nil meter probability)))
