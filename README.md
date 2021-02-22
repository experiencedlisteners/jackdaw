# Jackdaw

Jackdaw is a Common Lisp framework for defining discrete dynamic Bayesian networks with deterministic constraints in a way that involves writing very little code.

The formalities underlying this framework are described in chapter three of [Van der Weij (2020)](#vdweij2020).

Jackdaw was inspired by the [IDyOM modeling framework](https://github.com/mtpearce/idyom) by Marcus Pearce. 
The sequence model probability distributions included in jackdaw make use of IDyOM's implementation of the PPM algorithm.

## Usage example

The code sample below defines a fully functional implementation of a meter perception model presented in chapter five of [Van der Weij (2020)](#vdweij2020).
This model is based closely on the model described by [Van der Weij, Pearce, and Honing (2017)](#vdweij2017).

```common-lisp
(ql:require #:jackdaw)
(cl:in-package #:jackdaw)

(defmodel rhythm (dynamic-bayesian-network)
  (ioi-domain meter-domain)
  ((M ; meter
      (^m)
      (cpt ())
      (persist $^m meter-domain))
   (D ; downbeat distance
      (^d ^p m)
      (ppms (m))
      (chain (loop for ioi in ioi-domain
		   collect (cons (+ $^p ioi)
				 (ensure-list $^d)))
	     $^p))
   (P0 ; initial phase
       (^p0 m)
       (uniform ())
       (persist $^p0 (loop for p below (car $m) collect p)))
   (P ; phase
      (^p p0 m d)
      (uniform ())
      (recursive $^p (list (mod (car $d) (car $m)))
		 (list $p0)))
   (I ; inter-onset interval
      (d ^p ^i)
      (uniform ())
      (if (inactive? $d) (list +inactive+)
	  (list (- (car $d) $^p))))))
```

To instantiate, estimate, and query the model, we can use the REPL.

First, we need to instantiate the model.

```common-lisp
JACKDAW> (defparameter *model*
           (make-instance 'rhythm
                          :ioi-domain '(1 2 3 4)
                          :meter-domain '((8 4) (6 4) (6 8))
                          :p0-observer #'first
                          :m-observer (lambda (m) (list (second m) (third m)))
                          :i-observer (lambda (m) (if (listp m) (fourth m) m))))
````

Let's create some toy data from which we can estimate this instance.

First, we'll create a utility function for annotating sequences of IOIs with metrical information.

```common-lisp
(defun annotate (iois meter phase-0)
  "Utility function for annotating a list of IOIs with initial phase and meter."
  (loop for ioi in iois collect (cons phase-0 (append meter (list ioi)))))
```

Now, we can easily jot down some toy data and estimate the model.

```common-lisp
JACKDAW> (let ((data (list (annotate (list +inactive+ 4 2 2 4 1 1 1 1 4) '(8 4) 0)
                           (annotate (list +inactive+ 3 1 1 1 2 1 3 3) '(6 8) 0)
                           (annotate (list +inactive+ 2 1 1 2 2 1 1 2 2 2 2 1 1 4) '(6 4) 0))))
           (observe *model* 'i 'm 'p0)  ; configure I and M to be observed in *MODEL*
           (estimate *model* data)) ; estimate *MODEL* from the data
```

The following illustrates how the model instance can be queried on the REPL.

```common-lisp
JACKDAW> (hide *model* 'm 'p0)
NIL
JACKDAW> (probability *model* (list +inactive+ 1))
0.31908745
T
JACKDAW> (term:table 
          (state-probability-table
           (marginalize 
            (posterior-distribution
             *model*
             (generate *model* (list +inactive+ 4 2 2 4)))
            '(m))
           'm) :column-width 12)
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

# Documentation

An introductory tutorial to using jackdaw can be found [here](https://github.com/experiencedlisteners/jackdaw-tutorial).

# References

<a id="vdweij2017">van der Weij, B., Pearce, M. T., and Honing, H. (2017). A probabilistic model of meter perception: Simulating enculturation. *Frontiers in Psychology*. 8:824. doi: [10.3389/fpsyg.2017.00824](https://dx.doi.org/10.3389/fpsyg.2017.00824)
<a id="vdweij2017">van der Weij, B. (2020). *Experienced listeners: Modeling the influence of long-term musical exposure on rhythm perception.* (Doctoral dissertation, Universiteit van Amsterdam, Amsterdam) [PDF](https://hdl.handle.net/11245.1/dd3e25aa-6006-486e-afcf-c0692e0afacd)
