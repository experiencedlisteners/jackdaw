(ql:quickload "fare-csv")
(ql:quickload "jackdaw")
(ql:quickload "cl-ansi-term")

(defparameter *pitch-names-with-sharps* '(c c# d d# e f f# g g# a a# b ))
(defparameter *pitch-names-with-flats* '(c db d eb e f gb g ab a bb b))

(defun pitch-name (tonic pitch &key (octave-size 12) exclude-octave)
  (let* ((octave (floor (/ pitch octave-size)))
	 (degree (mod pitch octave-size))
	 (name
	   (if (<= tonic 7)
	       (elt *pitch-names-with-flats* degree)
	       (elt *pitch-names-with-sharps* degree))))
    (intern (format nil "~a~a" name (if exclude-octave "" octave)))))

(jd:defmodel tonality (jd:dynamic-bayesian-network)
  (&key order (octave 12) (pitch-alphabet (loop for p below 120 collect p)))
  ((tonic                                                ; the tonic (0-11)
	  (^tonic)                                       ; parents: conditional upon itself in the previous moment
	  (jd:uniform ())                                ; distribution: uniform
	  (jd:persist $^tonic                            ; constraint: generated only in the first moment and deterministic thereafter
		      (loop for tonic below octave collect tonic))
	  :observer #'first)
   (tonic-name ; tonic name for convenience
	       (Tonic)                             
	       (jd:uniform ())
	       (list
		(pitch-name $tonic $tonic
			    :octave-size octave
			    :exclude-octave t))
	       :observer #'first)
   (mode ; mode (major, minor)
	 ()                                              ; parents: none           
	 (jd:cpt ())                                     ; distribution: conditional probability table
	 (list 'major 'minor)                            ; constraint: possible values are 'major and 'minor
	 :observer #'second)                             ; observer: second
   (scale-degree ; scale degree
		 (^scale-degree mode)                    ; parents: itself in the previous moment, mode
		 (jd:ppms (mode) :order-bound order)     ; distribution: ppm model conditioned on mode
		 (jd:markov order $^scale-degree         ; constraint: markov model on scale degree alphabet = 0-11
			    (loop for deg below octave collect deg))
		 :observer #'third)                      ; observer: third
   (interval ; pitch interval
	     (^interval pitch ^pitch mode)
	     (jd:ppms (mode) :order-bound order)
	     (jd:markov order $^interval
			(if (jd:inactive? $^pitch) (list jd:+inactive+)
			    (list (cons $mode (- $pitch $^pitch)))))
	     :observer
	     (lambda (m)
	       (if (listp m) (third m) m)))
   (pitch ; chromatic pitch
	  (tonic scale-degree)                           ; parents: tonic and scale-degree
	  (jd:uniform ())                                ; distribution: uniform
	  (loop for pitch in pitch-alphabet              ; constraints: must be consistent with scale degree
		if (eq (mod (- pitch $tonic) 12)
		       (car $scale-degree))              
		  collect pitch)
	  :observer
	  (lambda (m)
	    (if (listp m) (third m) m)))
   (pitch-name ; pitch name for convenience
	       (tonic pitch)
	       (jd:uniform ())
	       (list
		(pitch-name $tonic $pitch
			    :octave-size octave))
	       :observer
	       (lambda (m)
		 (if (listp m) (third m) m)))))

(defun annotate-melody (keysig mode melody)
  (let ((scale (case mode
		 (0 'major)
		 (9 'minor)
		 (t (error "Unrecognized model ~a" mode))))
	(tonic (mod (* keysig 7) 12)))
    (loop for pitch in melody collect (list tonic scale pitch))))

(defun preprocess (rows)
  (let* ((dataset)
	 (rows (cdr rows))) ; skip header
    (dolist (row rows)
      (destructuring-bind (name keysig mode melody)
	  (mapcar #'read-from-string row)
	(declare (ignore name))
	(push (annotate-melody keysig mode melody) dataset)))
    (reverse dataset)))

(defun load-data (path)
  (let ((data (fare-csv:with-rfc4180-csv-syntax ()
		(fare-csv:read-csv-file path))))
    (preprocess data)))

(defmethod parameterize-model ((model tonality) path)
  (jd:hide model)
  (jd:observe model 'pitch 'mode 'tonic)
  (let ((data (load-data path)))
    (jd:estimate model data)))

;; Example usage 
;;
;; (defparameter *tonality* (make-musical-key-model))
;; (parameterize-model *tonality* "path/to/jackdaw-tutorial/materials/mtc-melodies.csv")
;; (jd:hide *tonality*)
;; (jd:observe *tonality* 'pitch)
;; (term:table 
;;  (jd:state-probability-table 
;;   (jd:posterior (jd:generate *tonality* '(7 5 4 0 2 7 0)))
;;   :variables '(tonic-name scale) :sort t)
;;  :column-width 15)
