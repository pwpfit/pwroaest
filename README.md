# `pwroaest` : Piecewise Region of Attraction Estimation

This toolbox offers sum-of-squares based algorithms for the region-of-attraction estimation for piecewise defined polynomial systems. Such systems are defined as ordinary differential equations `dx = f(x)` where
```
    f(x) = f1(x), if φ(x) <= 0
    f(x) = f2(x), else
```
and `f1, f2, φ` are polynomials in n variables `x`.

It is assumed throughout this work that `f1(0) = 0 if φ(0) <= 0` and `f2(0) = 0 if φ(0) >= 0`.


## Theoretical background

A set `X ⊂ ℝ^n` is an *estimate of the region of attraction of the equilibrium `x* = 0` for `f`* if any trajectory starting in `X` stays in `X` as well as converges to `x*` when time approaches infinity; i.e., `X` is both forward-invariant for `f` and attractive to `x*`.

### Definition: Lyapunov Functions
Let `f` be a system with `f(0) = 0`, `V` be a scalar field over `ℝ^n` and `X` a subset of `ℝ^n`; we will call `V` *Lyapunov function of `f` in `X`* if and only if
1. `V(0) = 0`
1. `V(x) > 0 for all x ∈ R^n - {0}`
1. `grad V(x)*f(x) < 0 for all x ∈ X - {0}`,

where `grad V(x)` denotes the gradient of `V` evaluated in `x`.

### Theorem: Piecewise Region of Attraction
Let `f1, f2, φ` define a piecewise system `f` as above, `X` be a subset of `ℝ^n`, and `γ > 0`; if one of the following hold,
- there is a scalar field `V` that is Lyapunov function of both `f1` and `f2` in `{x| V(x) <= γ}`, respectively, ("Local Lyapunov function"); and `X = {x| V(x) < γ}`
- there is a scalar field `V` that is Lyapunov function of both `f1` and `f2` in `{x| V(x) <= γ, φ(x) <= 0}` and `{x| V(x) <= γ, φ(x) >= 0}`, respectively, ("Common Lyapunov function"); and `X = {x| V(x) < γ}`
- there are scalar fields `V1, V2` such that `V1` is Lyapunov function of `f1` in `{x| V1(x) <= γ, φ(x) <= 0}` and `V2` is Lyapunov function of `f2` in `{x| V2(x) <= γ, φ(x) >= 0}`; `V1(x) = V2(x) for all x s.t. φ(x) = 0`; and `X = {x| V(x) < γ} ∩ {x| V(x) < γ}`

then `X` is an estimate of the region of attraction for `f`.
