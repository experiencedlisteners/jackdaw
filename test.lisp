(ql:quickload "jackdaw")

(jd::defdistribution meter
    (jd:cpt) (&key normalization-factor meter-cpt) (meter)
  (let ((meter-probability
	  (jd:probability meter-cpt (cons meter nil))))
    (pr:mul (car meter) (pr:div meter-probability normalization-factor))))

(jd::defestimator
    meter (data distribution) (meter) ()
    ((meter-cpt (jd:estimate (jd::make-cpt-distribution) data))
     (normalization-factor
      (apply #'pr:add
	     (loop
	       for meter in (jd:domain (meter-cpt distribution))
	       collect
	       (pr:mul (car meter)
		       (jd:probability (meter-cpt distribution) (cons meter nil))))))))

(jd:defmodel rhythm (jd:dynamic-bayesian-network)
  (ioi-domain meter-domain)
  ((M                    ; meter
      (^m)
      (jd:cpt ())   ; conditional probability table
      (jd:persist $^m meter-domain))
   (D                    ; downbeat distance
      (^d ^p m)
      (jd:ppms (m)) ; set of PPM sequence models
      (jd:chain (loop for ioi in ioi-domain
			   collect (cons (+ $^p ioi)
					 (jd:ensure-list $^d)))
		     $^p))
   (P0                   ; initial phase (or pickup interval)
       (^p0 m)
       (jd:uniform ())
       (jd:persist $^p0 (loop for p below (car $m) collect p)))
   (P                    ; phase
      (^p p0 m d)
      (jd:uniform ())
      (jd:recursive $^p (list (mod (car $d) (car $m)))
			 (list $p0)))
   (I                    ; inter-onset interval
      (d ^p ^i)
      (jd:uniform ())
      (if (jd:inactive? $d) (list jd:+inactive+)
	  (list (- (car $d) $^p))))))
