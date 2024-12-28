@doc raw"""
    LinearDispersionRelation(ref_height)

A struct representing a linear dispersion relation ``\omega(k)`` of an equation. The reference
water height `h0` is given by `ref_height`.
"""
struct LinearDispersionRelation{RealT <: Real}
    ref_height::RealT
end

function Base.show(io::IO, disp_rel::LinearDispersionRelation)
    print(io, "LinearDispersionRelation(h0 = ", disp_rel.ref_height, ")")
end

Base.broadcastable(disp_rel::LinearDispersionRelation) = Ref(disp_rel)

@doc raw"""
    wave_speed(disp_rel, equations, k; normalize = false)

Compute the wave speed ``c`` for a given wavenumber ``k`` using the linear dispersion relation ``disp_rel``
of the `equations`.
The wave speed is given by ``c = \omega(k) / k``. If `normalize` is `true`, the wave speed is normalized
by the shallow water wave speed ``\sqrt{g h0}``, where `g` is the gravity constant and `h0` is the reference
water height of the dispersion relation.
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

# TODO: Factor of 5/2 in front?
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
# - Didier Clamond Denys Dutykh, Dimitrios Mitsotakis (2017)
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
