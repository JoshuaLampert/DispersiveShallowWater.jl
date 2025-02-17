@doc raw"""
    LinearDispersionRelation(ref_height)

A struct representing a linear dispersion relation ``\omega(k)`` of an equation. The reference
water height `h0` is given by `ref_height`. A dispersion relation can be called as
`disp_rel(equations, k)` to compute the wave frequency ``\omega(k)`` for a given wavenumber `k`
and a set of equations.

See also [`wave_speed`](@ref) for computing the wave speed ``c = \omega(k) / k`` given a linear
dispersion relation.
"""
struct LinearDispersionRelation{RealT <: Real}
    ref_height::RealT
end

function Base.show(io::IO, disp_rel::LinearDispersionRelation)
    print(io, "LinearDispersionRelation(h0 = ", disp_rel.ref_height, ")")
end

Base.broadcastable(disp_rel::LinearDispersionRelation) = (disp_rel,)

@doc raw"""
    wave_speed(disp_rel, equations, k; normalize = false)

Compute the wave speed ``c`` for a given wavenumber ``k`` using the
[`LinearDispersionRelation`](@ref) `disp_rel` of the `equations`.
The wave speed is given by ``c = \omega(k) / k``. If `normalize` is `true`, the wave speed is normalized
by the shallow water wave speed ``\sqrt{g h_0}``, where ``g`` is the `gravity_constant` of the `equations`
and ``h_0`` is the `ref_height` of the dispersion relation `disp_rel`.

See also [`LinearDispersionRelation`](@ref).
"""
function wave_speed(disp_rel::LinearDispersionRelation, equations, k;
                    normalize = false)
    omega = disp_rel(equations, k)
    c = omega / k
    if normalize
        c /= sqrt(gravity_constant(equations) * disp_rel.ref_height)
    end
    return c
end

@doc raw"""
    EulerEquations1D(; gravity_constant, eta0 = 0.0)

A struct representing the 1D Euler equations with a given gravity constant and a still-water surface
`eta0`.

!!! note
    In DispersiveShallowWater.jl, the Euler equations are *only* used for computing the full linear dispersion
    relation
    ```math
    \omega(k) = \sqrt{g k \tanh(h_0 k)}.
    ```
    They cannot be solved as a system of equations.
"""
struct EulerEquations1D{RealT <: Real} <: AbstractEquations{1, 0}
    gravity::RealT
    eta0::RealT
end

function EulerEquations1D(; gravity_constant, eta0 = 0.0)
    return EulerEquations1D(gravity_constant, eta0)
end

function (disp_rel::LinearDispersionRelation)(equations::EulerEquations1D, k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    return sqrt(g * k * tanh(h0 * k))
end

function (disp_rel::LinearDispersionRelation)(equations::BBMEquation1D, k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    return sqrt(g * h0) * k / (1 + 1 / 6 * (h0 * k)^2)
end

# See
# - Joshua Lampert, Hendrik Ranocha (2024)
#  Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#  [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
function (disp_rel::LinearDispersionRelation)(equations::BBMBBMEquations1D, k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    return sqrt(g * h0) * k / (1 + 1 / 6 * (h0 * k)^2)
end

# See
# - Magnus Svärd, Henrik Kalisch (2023)
#   A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
#   [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)
function (disp_rel::LinearDispersionRelation)(equations::SvärdKalischEquations1D, k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    c0 = sqrt(g * h0)
    alpha = equations.alpha * c0 * h0^2
    beta = equations.beta * h0^3
    gamma = equations.gamma * c0 * h0^3
    a = (1 + beta / h0 * k^2)
    b = (-alpha - beta * alpha / h0 * k^2 - gamma / h0) * k^3
    c = -g * h0 * k^2 + gamma * alpha / h0 * k^6
    return (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
end

# See, e.g., eq. (45) (α = 1) in
# - Maria Kazolea, Andrea Filippini, Mario Ricchiuto (2023)
#   Low dispersion finite volume/element discretization of the enhanced Green-Naghdi equations for
#   wave propagation, breaking and runup on unstructured meshes
#   [DOI: 10.1016/j.ocemod.2022.102157](https://doi.org/10.1016/j.ocemod.2022.102157)
#
# or eq. (53) (β = 0) in
# - Didier Clamond, Denys Dutykh, Dimitrios Mitsotakis (2017)
#   Conservative modified Serre–Green–Naghdi equations with improved dispersion characteristics
#   [DOI: 10.1016/j.cnsns.2016.10.009](https://doi.org/10.1016/j.cnsns.2016.10.009)
function (disp_rel::LinearDispersionRelation)(equations::SerreGreenNaghdiEquations1D, k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    return sqrt(g * h0) * k / sqrt(1.0 + (k * h0)^2 / 3)
end

# TODO: Make this for general lambda (how to understand eq. (19) in Favrie and Gavrilyuk?)
function (disp_rel::LinearDispersionRelation)(equations::HyperbolicSerreGreenNaghdiEquations1D,
                                              k)
    h0 = disp_rel.ref_height
    g = gravity_constant(equations)
    lambda = equations.lambda
    return sqrt(g * h0) * k / sqrt(1.0 + (k * h0)^2 / 3)
end
