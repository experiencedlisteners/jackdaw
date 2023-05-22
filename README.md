# Jackdaw

Jackdaw is a Common Lisp framework for defining discrete dynamic Bayesian networks with deterministic constraints in a way that involves writing very little code.

The formalities underlying this framework are described in chapter three of [Van der Weij (2020)](#vdweij2020).

Jackdaw was inspired by the [IDyOM modeling framework](https://github.com/mtpearce/idyom) by Marcus Pearce. 
The sequence model probability distributions included in jackdaw make use of IDyOM's implementation of the PPM algorithm.

## Usage example

The code sample below defines a fully functional implementation of a meter perception model presented in chapter five of [Van der Weij (2020)](#vdweij2020).
This model is based closely on the model described by [Van der Weij, Pearce, and Honing (2017)](#vdweij2017).

```common-lisp
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
```

The above first uses `DEFDISTRIBUTION` to defin a custom categorical probability distribution (which uses the built-in `JACKDAW:CPT`, an implementation of conditional probability tables).
This custom distribution incorporates a factor that compensates for the fact that, in the world of discrete symbolic rhythms, metrical interpretations with longer periods have to spread out probability mass over more possible initial phases. This factor is calculated in the custom estimator defined below with `DEFESTIMATOR`.

Next, `DEFMODEL` defines a dynamic Bayesian network graph with five variables (`M`, `D`, `P0`, `P`, and `I`), their graphical dependency relations, their probability distributions, and their congruency constraints.
For a detailed description of this model, see chapter five of [Van der Weij (2020)](#vdweij2020).
For the purpose of this example it is useful to know that this model is generates sequences of inter-onset intervals.
The first inter-onset interval in such sequences must always be the constant `+INACTIVE+`, defined in the `jackdaw` package, since this moment is used to generate the initial phase (or pickup interval).

We can use the REPL to instantiate, estimate, and query the model.

Instantiating the model could be done as follows.

```common-lisp
CL-USER> (defparameter
          *model*
          (make-instance 'rhythm
                         :ioi-domain '(1 2 3 4)
                         :meter-domain '((8 4) (6 4) (6 8))
                         :p0-observer #'first
                         :m-observer (lambda (m) (list (second m) (third m)))
                         :i-observer (lambda (m) (if (listp m) (fourth m) m))))
````

We have to provide values for the model's two parameters: `IOI-DOMAIN` and `METER-DOMAIN`.
These should be interpreted as lists of inter-onset intervals and metrical interpretations that the model can generate.

In order to estimate the model, let's create some toy data.

First, we'll create a utility function for annotating sequences of IOIs with metrical information.

```common-lisp
(defun rhythm (iois &optional meter phase-0)
  "Utility function for annotating a list of IOIs with initial phase and meter."
  (flet ((annotate (ioi)
	   (cond ((and (null meter) (null phase-0))
		  ioi)
		 ((null phase-0)
		  (error "Providing only a meter and no initial phase is not allowed"))
		 (t
		  (cons phase-0 (append meter (list ioi)))))))
    (mapcar #'annotate (cons jd:+inactive+ iois))))
```

Now, we can easily jot down some data and estimate the model.

```common-lisp
CL-USER> (let ((data (list (rhythm '(4 2 2 4 1 1 1 1 4)           '(8 4) 0)
                           (rhythm '(3 1 1 1 2 1 3 3)             '(6 8) 0)
                           (rhythm '(2 1 1 2 2 1 1 2 2 2 2 1 1 4) '(6 4) 0))))
           (jd:hide *model*) ; hide the entire model
           (jd:observe *model* 'i 'm 'p0)  ; configure I and M to be observed in *MODEL*
           (jd:estimate *model* data)) ; estimate the model from the data
```

Above, we first used `OBSERVE` to make the variables `I`, `M`, and `P0` of the model observable.
Providing values for these variables is sufficient to make the model fully observable.
Then we used `ESTIMATE` to estimate the model from the data

The following illustrates how the model instance can be queried on the REPL.

```common-lisp
CL-USER> (jd:hide *model*) ; hide all variables
NIL
CL-USER> (jd:observe *model* 'i) ; observe inter-onset interval
NIL
CL-USER> (jd:probability *model* (rhythm '(1)))
0.31908745
T
CL-USER> (ql:quickload "cl-ansi-term") ; for nice output formatting
...
CL-USER> (term:table 
          (jd:state-probability-table
           (jd:marginalize 
            (jd:posterior
             (jd:generate *model* (rhythm '(4 2 2 4))))
           '(m))) :column-width 12)
+-----------+-----------+
|M          |PROBABILITY|
+-----------+-----------+
|(8 4)      |0.7039645  |
+-----------+-----------+
|(6 4)      |0.24176745 |
+-----------+-----------+
|(6 8)      |0.0542681  |
+-----------+-----------+
NIL
```

Above, we first used `HIDE` to make `M` and `P0`, the two variables that describe the metrical interpretation of a rhythm, hidden.

Then, we evaluated the *model evidence* for an (otherwise uninteresting) rhythmic pattern with `PROBABILITY`.

Finally, we generated model states congruent with a simple rhythmic pattern, defined by the list of inter-onset intervals `(4 2 2 4)` with `GENERATE`.
We marginalized them to a states containing just the variable `M` with `MARGINALIZE`, calculated posterior probabilities with `POSTERIOR`.
The result of these operations is a list of states, which we converted to a probability table with `STATE-PROBABILITY-TABLE`, and displayed using `TABLE` from the `cl-ansi-term` library.

# Documentation

An introductory tutorial to using jackdaw can be found [here](https://github.com/experiencedlisteners/jackdaw-tutorial).

# References

<a id="vdweij2017">van der Weij, B., Pearce, M. T., and Honing, H. (2017). A probabilistic model of meter perception: Simulating enculturation. *Frontiers in Psychology*. 8:824. doi: [10.3389/fpsyg.2017.00824](https://dx.doi.org/10.3389/fpsyg.2017.00824)

<a id="vdweij2020">van der Weij, B. (2020). *Experienced listeners: Modeling the influence of long-term musical exposure on rhythm perception.* (Doctoral dissertation, Universiteit van Amsterdam, Amsterdam) [PDF](https://hdl.handle.net/11245.1/dd3e25aa-6006-486e-afcf-c0692e0afacd)
