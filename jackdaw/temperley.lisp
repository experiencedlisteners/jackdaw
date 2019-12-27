(cl:in-package #:jackdaw)

;;;;;;;;;;;;;;;;;;;; Probability distributions ;;;;;;;;;;;;;;;;;;;;

(defclass tactus-interval (distribution)
  ((t0 :initarg :t0 :accessor t0)))
(defclass beat-deviation (distribution)
  ((subdivision :initarg :subdivision :reader subdivision)
   (phase :initarg :phase :reader beat-phase :initform 1)
   (p :initarg :p :reader p)))
(defclass bar-phase (distribution)
  ((duple :initarg :duple :reader duple :initform .5)
   (triple :initarg :triple :reader triple :initform '(.3 .3))))
(defclass tactus-phase (distribution)
  ((zero :initarg :zero :reader zero
	 :documentation
	 "Probability of a phase of zero. Must be > 0 and <= 1.")))
(defclass note (distribution)
  ((p :initarg :p :reader p :initform '(.5 .5 .5 .5))))
(defclass temperley-ioi (note) ())

;; Serializers

(defwriter tactus-interval (m) (t0 m))
(defreader tactus-interval (m t0)
  (setf (slot-value m 't0) t0))
(defreader beat-deviation (m d)
  (setf (slot-value m 'subdivision) (first d))
  (setf (slot-value m 'phase) (second d))
  (setf (slot-value m 'p) (third d)))
(defwriter beat-deviation (m) (list (subdivision m) (phase m) (p m)))
(defreader beat-deviation (m d)
  (setf (slot-value m 'subdivision) (first d))
  (setf (slot-value m 'phase) (second d))
  (setf (slot-value m 'p) (third d)))
(defwriter bar-phase (m) (list (duple m) (triple m)))
(defreader bar-phase (m d)
  (setf (slot-value m 'duple) (first d))
  (setf (slot-value m 'triple) (second d)))
(defwriter tactus-phase (m) (zero m))
(defreader tactus-phase (m zero) (setf (slot-value m 'zero) zero))
(defwriter note (m) (p m))
(defreader note (m p) (setf (slot-value m 'p) p))

(defmethod probability ((d tactus-interval) arguments symbol)
  (let ((previous (car arguments)))
    (pr:in (if (eq previous +inactive+)
	       (elt (t0 d) (- symbol 9)) ;(min-tactus (model d)))
	       (exp (- (expt (* 0.5 (- symbol previous )) 2)))))))

(defmethod probability ((d note) arguments symbol)
  (let ((p (elt (p d) (car arguments))))
    (pr:in (if symbol p (- 1 p)))))

(defmethod probability ((d bar-phase) arguments symbol)
  "Calculate the probability of a tactus beat phase with respect to the bar.
Arguments is a list of length one containing the number of tactus beats per bar."
  (let ((grouping (car arguments)))
    (pr:in
     (case grouping
       (2 (elt (append (duple d) (list (- 1 (apply #'+ (duple d))))) symbol))
       (3 (elt (append (triple d) (list (- 1 (apply #'+ (triple d))))) symbol))))))

(defmethod probability ((d beat-deviation) arguments symbol)
  (if (eq symbol '*) (pr:in 1)
      (let* ((period (car arguments))
	     (center (beat-location period (subdivision d) (beat-phase d)))
	     (deviation (abs (- center symbol))))
	(pr:in (elt (p d) deviation)))))

(defmethod probability ((d tactus-phase) arguments symbol)
  (let ((interval (car arguments))
	(p-zero (zero d)))
    (pr:in
     (if (eq symbol 0) p-zero
	 (/ (- 1 p-zero) (1- interval))))))
	
(defun beat-location (period subdivision phase)
  (floor (/ (* phase period) subdivision)))

(defun beat-locations (max-beat-dev s period phase &optional (previous 0))
  (let ((center (beat-location period s phase)))
    (loop for p
       from (max (1+ previous) (- center max-beat-dev))
       to (min (- period (- s phase)) (+ center max-beat-dev))
       collect p)))

;;;;;;;;;;;;;;;;;;;; Congruency constraints ;;;;;;;;;;;;;;;;;;;;

(defmodel temperley (generative-model)
  ((observed :initform '(:n))
   (tacti :initarg :tacti :reader tacti
	  :initform (loop for $t from 9 below 23 collect $t))
   (max-beat-dev :initarg :max-beat-dev :reader max-beat-dev :initform 0));3))
  ((U (^u) (one-shot '(2 3)))
   (L (^l) (one-shot '(2 3)))
   (T (^t ^tph)
       (recursive (if (eq (1+ $^tph) $^t) (tacti model) (list $^t))
		  (tacti model)))
   (BPH (^bph tph u)
	(recursive (if (eq $tph 0)
		       (list (mod (1+ $^bph) $u))
		       (list $^bph))
		   (loop for bph below $u collect bph)))
   (TPH (^tph ^t t)
	(recursive (list (mod (1+ $^tph) $^t))
		   (loop for tph below $t collect tph)))
   (DB (^db tph t l)
       (recursive (cond ((eq $l 3) +singleton+)
			((eq $tph 0) (beat-locations (max-beat-dev model) 2 $t 1))
			(t (list $^db)))
		  (if (eq $l 3) +singleton+ (beat-locations (max-beat-dev model) 2 $t 1))))
   (TB1 (^tb1 tph t l)
	(recursive (cond ((eq $l 2) +singleton+)
			 ((eq $tph 0) (beat-locations (max-beat-dev model) 3 $t 1))
			 (t (list $^tb1)))
		   (if (eq $l 2) +singleton+ (beat-locations (max-beat-dev model) 3 $t 1))))
   (TB2 (^tb2 tb1 tph t l)
	(recursive (cond ((eq $l 2) +singleton+)
			 ((eq $tph 0) (beat-locations (max-beat-dev model) 3 $t 2 $tb1))
			 (t (list $^tb2)))
		   (if (eq $l 2)
		       +singleton+
		       (beat-locations (max-beat-dev model) 3 $t 2 $tb1))))
   (BS (db tb1 tb2 tph bph)
       (normal (list (cond ((and (eq $tph 0) (eq $bph 0)) 3)
			   ((eq $tph 0) 2)
			   ((member $tph (list $db $tb1 $tb2)) 1)
			   (t 0)))))
   (N (^n) (recursive (list t nil) (list t))))
  ((U nil (bernouilli :symbols '(3 2)))
   (L nil (bernouilli :symbols '(3 2)))
   (T (^t) (tactus-interval))
   (BPH (u) (bar-phase))
   (TPH (t) (tactus-phase))
   (DB (t) (beat-deviation :subdivision 2))
   (TB1 (t) (beat-deviation :subdivision 3))
   (TB2 (t) (beat-deviation :subdivision 3 :phase 2))
   (N (bs) (note))))

(defmethod initialize-instance ((m temperley)
				&key (u .24) (l .22)
				  (t0 '(.1 .2 .3 .23 .13 .03 .006 .002 .001
					.0006 .0002 .0001 .00005 .00005))
				  (duple-ph0 .65)
				  (triple-ph0 .33)
				  (triple-ph1 .667)
				  (t-ph0 .6)
				  (deviation '(.32 .24 .08 .02))
				  (note '(.01 .38 .74 .95)))
  "U is the probability of the upper level being triple, L of the lower level 
begin triple, T0 a list of probabilities of initial tactus intervals (from
minimum to maximum), DUPLE-PH0 the probability of a BPH of zero given a duple 
upper level, TRIPLE-PH0 and TRIPLE-PH1 are the probabilities of BPH being 1 or 2 
given a triple upper level. T-PH0 is the probability of the tactus phase being zero. 
NOTE is a list of probabilities of note onsets at metrical levels from 0 (not a beat) to
the maximum beat salience. DEVIATION is a list of probabilities that a note deviates from
a beat by N pips, where N is the position in the list."
  (setf (slot-value (distribution m 'U) 'p) u
	(slot-value (distribution m 'L) 'p) l
	(slot-value (distribution m 'T) 't0) t0
	(slot-value (distribution m 'BPH) 'duple) (list duple-ph0)
	(slot-value (distribution m 'BPH) 'triple) (list triple-ph0 triple-ph1)
	(slot-value (distribution m 'TPH) 'zero) t-ph0
	(slot-value (distribution m 'DB) 'p) deviation
	(slot-value (distribution m 'TB1) 'p) deviation
	(slot-value (distribution m 'TB2) 'p) deviation
	(slot-value (distribution m 'N) 'p) note))
						
(defmodel event-based-symbolic-temperley (generative-model)
  ((observed :initform '(:n))
   (tacti :reader tacti)
   (ioi-domain :initarg :ioi-domain :reader ioi-domain))
  ((U (^u) (one-shot '(2 3)) :inputs (u) :posterior-constraint (eq $u <u))
   (L (^l) (one-shot '(2 3)) :inputs (l) :posterior-constraint (eq $l <l))
   (T (^t) (one-shot (tacti model))
      :inputs (t) :posterior-constraint (eq $t <t)) ; changed to one-shot
   (BPH (^bph u) (one-shot (loop for bph below $u collect bph))
	 :inputs (bph) :posterior-constraint (eq $bph <bph))
   (TPH (^tph t) (one-shot (loop for tph below $t collect tph))
	      :inputs (tph) :posterior-constraint (eq $tph <tph))
   (PH (^ph ioi u t bph tph)
       (recursive (deterministic (mod (+ $^ph $ioi) (* $u $t)))
		  (deterministic (+ $tph (* $bph $t)))))
   ;;:inputs (ioi)
   ;;:posterior-constraint
   ;;(recursive t (eq $ph <ioi)))
   (IOI (^ph)
	(chain (ioi-domain model) (^ph))
	:inputs (ioi)
	:posterior-constraint
	(chain-posterior (eq <ioi $ioi) (^ph))))
  ((U () (bernouilli :symbols '(3 2)))
   (L () (bernouilli :symbols '(3 2)))
   (T () (categorical))
   (BPH (u) (bar-phase))
   (TPH (t) (tactus-phase))
   (IOI (^ph t u l) (temperley-ioi))))

(defmethod initialize-instance :after ((m event-based-symbolic-temperley)
				&key u l t0 t-ph0 duple-ph0 triple-ph0 triple-ph1 note)
  "T0 must be a list of pairs containing tactus and corresponding probability."
  (setf (slot-value (distribution m 'U) 'p) u
	(slot-value (distribution m 'L) 'p) l
	(slot-value (distribution m 'BPH) 'duple) (list duple-ph0)
	(slot-value (distribution m 'BPH) 'triple) (list triple-ph0 triple-ph1)
	(slot-value (distribution m 'TPH) 'zero) t-ph0
	(slot-value (distribution m 'IOI) 'p) note
	(slot-value m 'tacti) (mapcar #'car t0))
  (hide m 'U 'L 'T 'BPH 'TPH)
  (loop for (interval p) in t0 do
       (set-param (distribution m 'T) nil interval p)))

(defun beat-salience ($ph $u $l $t)
  (cond ((eq (mod $ph (* $u $t)) 0) 3)
	((eq (mod $ph $t) 0) 2)
	((eq (mod $ph (/ $t $l)) 0) 1)
	(t 0)))

(defmethod probability ((d temperley-ioi) arguments ioi)
  (let (($^ph (first arguments)) ($t (second arguments))
	($u (third arguments)) ($l (fourth arguments)))
    ;; Product of probability of onset at SYMBOL
    ;; and probabilities of no onsets at positions from PREVIOUS + 1 to
    ;; SYMBOL - 1.
    (pr:in (* (apply #'* (loop for interval from 1 below ioi collect
			      (- 1 (elt (p d) (beat-salience (+ $^ph interval) $u $l $t)))))
	      (elt (p d) (beat-salience (+ $^ph ioi) $u $l $t))))))

(defmethod moment ((m temperley) moment congruent-states)
  ;; Todo implement more efficient state generation
  (call-next-method))

(defmethod generate-states ((m temperley) vertices previous-state moment &optional root?)
  ;; Todo implement more efficient state generation
  (call-next-method))

(defwriter event-based-symbolic-temperley (m)
  (list (read-from-string (with-output-to-string (s) (call-next-method m s)))
	(tacti m) (ioi-domain m)))
(defreader event-based-symbolic-temperley (m data)
  (let ((model (first data))
	(tacti (second data))
	(ioi-domain (third data)))
    (with-input-from-string (s (write-to-string model))
      (call-next-method m s)
      (setf (slot-value m 'ioi-domain) ioi-domain)
      (setf (slot-value m 'tacti) tacti))))
