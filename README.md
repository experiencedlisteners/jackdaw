# Jackdaw

Jackdaw is a Common Lisp framework for defining discrete dynamic Bayesian networks with deterministic constraints in a way that involves writing very little code.

The formalities underlying this framework are described in chapter three of my [PhD thesis](https://dare.uva.nl/search?identifier=dd3e25aa-6006-486e-afcf-c0692e0afacd).

Below is an example of a full definition in jackdaw of a meter perception model presented in chapter five of my PhD thesis.
This model is based closely on the model described by [Van der Weij, Pearce, and Honing (2017)](#vdweij2017).

```common-lisp
(in-package #:jackdaw)

(defmethod meter-domain ((m enculturation))
  (loop for m in (meter-parameters m) collect (car (car m))))

(defmodel enculturation ()
  (training? ioi-domain meter-parameters) 
  ((M ; meter
      (^m)
      (categorical
       () :parameters (meter-parameters enculturation))
      (persistent (meter-domain model))
      :key (lambda (m) (list (getf m 'period) (getf m 'pulses))))
   (D ; downbeat distance
      (^d ^p m)
      (accumulator-model (m)
			 :training? (training? enculturation))
      (recursive (loop for ioi in (ioi-domain model)
		       collect (cons (+ $^p ioi) $^d))
		 (deterministic '(*)))
      :formatter (lambda (d) (car d)))
   (P ; phase
      (^p m d)
      (uniform ())
      (recursive (list (mod (car $d) (car $m)))
		 (loop for phase below (car $m)
		       collect phase)))
   (IOI ; inter-onset interval
        (d ^p ^ioi)
        (uniform ())
	    (if (eq (car $d) '*) '(*)
	        (list (- (car $d) $^p))))))
```

# Documentation

An introductory tutorial can be found [here](https://github.com/experiencedlisteners/jackdaw-tutorial).

# Examples

- model examples
- cli examples

# Related links

# Models

Examples of jackdaw models can be found in [this repository](https://github.com/experiencedlisteners/jackdaw-models).
- jackdaw cli
- data notebooks?

# Acknowledgments

The creation of Jackdaw was inspired by the [IDyOM framework](https://github.com/mtpearce/idyom) by Marcus Pearce. 
The sequence model probability distributions included in jackdaw make use of IDyOM's implementation of the PPM algorithm.

# References

<a id="vdweij2017">van der Weij, B., Pearce, M. T., and Honing, H. (2017). A probabilistic model of meter perception: Simulating enculturation. *Frontiers in Psychology*. 8:824. doi: [10.3389/fpsyg.2017.00824](https://dx.doi.org/10.3389/fpsyg.2017.00824)
