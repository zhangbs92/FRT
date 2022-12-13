# Functional regression based test for common and rare variants
User friendly software that apply functional regression test to genetic data

This is the functional regression based test (FRT) repo for genetic variant association test. Rather than seeing a variant as random variable, FRT see them as a stochastic process. Therefore, the genetic effect is a continuous function that changing with the relative position of a variant on a sequence. 

Currently this method can be applied to analyze:
a. population quantitative traits
b. family quantitative traits
c. population binary traits
d. family binary traits
e. population survival traits
f. related sample survival traits
g. population related survival traits
h. longitudinal quantitative traits
i. longitudinal binary traits

All 9 methods were implement in R by different codes, this repo is to unify everything in the same framework, hoping to extent the method into other field.
