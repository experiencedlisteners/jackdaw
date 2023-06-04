

(defun three-of-a-kind (roll)
  

(jd:defmodel yatzee (bayesian-network) ()
  ((dice1 () (jd:uniform '(1 2 3 4 5 6)))
   (dice2 () (jd:uniform '(1 2 3 4 5 6)))
   (dice3 () (jd:uniform '(1 2 3 4 5 6)))
   (dice4 () (jd:uniform '(1 2 3 4 5 6)))
   (dice5 () (jd:uniform '(1 2 3 4 5 6)))
   (roll (dice1 dice2 dice3 dice3 dice4 dice5) (jd:uniform)
	 (list dice1 dice2 dice3 dice4 dice5))
   (category (roll)
	     (loop 
