
**Simple distribution**

# Beta

Anything that are intrinsically bounded by 0 and 1. Emphasizing "intrinsical" is to set it apart from `logit_normal` distribution where it is also bounded by 0 and 1 due to the logi function. But it is more like a continuous random variable whereas a bound has been added.

- You can use `Beta (s,s)` to represent distribution peaked at 0.5, `s` determine how informative this disitribution will be

- You can use `Beta (0.5,0.5)` to represent bimodal distribution

- When `a <= 1`, it is left skewed, when `b <= 1`, it is right skewed.

# Gamma

It is very flexible, anything that is positive can be modeled by Gamma, you just need to change the parameter to flexibly adjust the variance and shape. `chi2` distribution is a special case for Gamma. 

- You can use `Gamma (m,1)` to make sure it is peaked at `m`

- You can use `Gamma (s,s)` to make sure it is peaked at 1, and `s` determine how much regularisation you'd like to add on that when deviating from 1.

# LogNorm

It is to measure any positive values, without clear overdispersion.

# Normal distribution

Honestly the default, like to model the coefficent of the regression.

# Student T distribution

Normal distribution with heavier tails

# Poisson

Default for counts data, the intuition is during certain time interval, how many buses arriving.

# Exponential

Related to the Poisson, the intuition is for the same situation in Poisson, the Exponential can measure the waiting time for the buses or the interval to be more general. It can be used for regularization as the probability goes down when the randome variable increases.

# Negative Bionomial

Poisson with overdispersion for counting data like RNA count, the intuition is, when achieving certain success, how many failures have been experienced.

# Binomial

Very intuitive, total number of trials each has a probability to succeed.

# Bernoulli

Very intuitive, toss the coin, a special case for binomial

# Categorical

Very intuitive, toss the dice, a special case for multinomial

# Multinomial

Very intuitive, toss n number of dice, the count of each side 

# Dirichlet

The multi-variate version of Beta, parameterizing Multinomial and categorical distribution

# Wishard

The multi-variate version of Gamme distribution.

# HyperGeometric

Very strong intuition, when having certain amount of fake product, and try to sample certain amount of product from total product.

**Compound distribution**

# Spike-and-Slab (or mixutre in general)

A mixed distribution of a spike and a slab, useful for variable selection, a spike is a distribution that puts all mass into a center, whereas slab is a distribution that gives mass to all values.

# Zero-inflated
Add inflated zero component to like Poisson or Negative Binomial

# Affine transformation

The loc and scale parameter




