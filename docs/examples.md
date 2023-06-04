# Examples

## Probability distributions (both distributions and jackdaw models)

- probability
- probabilities

- defreader
- defwriter
- probabilities


## Distributions

The `DISTRIBUTIONS` package provides tools for defining custom probability distributions as well as a number of pre-defined probability distributions.

The distinction between probability distributions and jackdaw models is in a way artificial, since jackdaw models can also be treated as (typically more complicated) probability distributions.

To define a new distribution, you have to define a class identifying the distribution and a method of calculating probabilities from parameters.

For convenience, you can use `DEFDISTRIBUTION`.

For defining the a conditional probability distribution, you have two options.

1. provide an implementation of the `PROBABILITY` method for your distribution
2. provide an implementaiton of the `PROBABILITIES` method for your distribution

Jackdaw only calls `PROBABILITIES`. The fallback behaviour of `PROBABILITIES` is to call `PROBABILITY` for each parameter in the domain, but you are free to provide your own implementation of `PROBABILITIES`.

API

methods

- estimate
- probability
- probabilities

macros

### DEFESTIMATOR


### DEFDISTRIBUTION



## Jackdaw

### Model definition

### Inference
