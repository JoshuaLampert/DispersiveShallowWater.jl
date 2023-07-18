@doc raw"""
    ShallowWaterEquations1D(gravity, eta0)

Basic non-dispersive shallow water equations (SWE) in one spatial dimension with spatially varying bathymetry. The equations are given by
```math
\begin{aligned}
  \frac{\partial\eta}{\partial t} + \frac{\partial}{\partial x}((\eta + D)v) &= 0,\\
  \frac{\partial v}{\partial t} + g\frac{\partial\eta}{\partial x} + \frac{\partial}{\partial x}\left(\frac{1}{2}v^2\right) &= 0.
\end{aligned}
```
The unknown quantities of the shallow equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography (bathymetry) ``b = -D > 0``. The water height above the bathymetry is therefore given by
``h = \eta + D``.

References for the SWE are many but a good introduction is available in Chapter 13 of the book:
- Randall J. LeVeque (2002)
  Finite Volume Methods for Hyperbolic Problems
  [DOI: 10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253)

"""
struct ShallowWaterEquations1D{RealT <: Real} <: AbstractShallowWaterEquations{1, 3}
    gravity::RealT # gravitational constant
    eta0::RealT    # constant "lake-at-rest" total water height
end

function ShallowWaterEquations1D(; gravity_constant, eta0 = 0.0)
    ShallowWaterEquations1D(gravity_constant, eta0)
end

varnames(::ShallowWaterEquations1D) = ("eta", "v", "D")

# TODO: Initial condition should not get a `mesh`
# TODO: Find good convergence test
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterEquations1D, mesh)


"""
function initial_condition_convergence_test(x,
                                            t,
                                            equations::ShallowWaterEquations1D,
                                            mesh)
    eta = 1.0
    v = 0.0
    D = 0.0
    return SVector(eta, v, D)
end

# TODO: Initial condition should not get a `mesh`
"""
    initial_condition_radial_dambreak(x, t, equations::ShallowWaterEquations1D, mesh)

An initial condition with a jump in the waterheight, still water and a Gaussian
bump as bathymetry.
"""
function initial_condition_radial_dambreak(x, t, equations::ShallowWaterEquations1D, mesh)
    eta = abs(x) < 0.5 ? 1.1 : 0.9
    v = 0.0
    D = -0.5 * exp(-5.0 * x^2)
    return SVector(eta, v, D)
end

function create_cache(mesh,
                      equations::ShallowWaterEquations1D,
                      solver::Solver,
                      initial_condition,
                      RealT,
                      uEltype)
    tmp1 = Array{RealT}(undef, nnodes(mesh))
    return (tmp1 = tmp1,)
end

# Split form discretization that conserves the mass (for eta and hv) and the energy for periodic boundary conditions, see
# - Niklas Wintermeyer, Andrew R. Winters , Gregor J. Gassner , David A. Kopriva (2016)
#   An Entropy Stable Nodal Discontinuous Galerkin Method for the Two
#   Dimensional Shallow Water Equations on Unstructured Curvilinear Meshes with
#   Discontinuous Bathymetry
#   [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
# - Hendrik Ranocha (2016)
#   Shallow water equations: split-form, entropy stable, well-balanced, and positivity
#   preserving numerical methods
#   [DOI: 10.1007/s13137-016-0089-9](https://doi.org/10.1007/s13137-016-0089-9)
# Here, in primitive variables instead of in conservative.
function rhs!(du_ode, u_ode, t, mesh, equations::ShallowWaterEquations1D,
              initial_condition,
              ::BoundaryConditionPeriodic, solver, cache)
    @unpack tmp1 = cache

    u = wrap_array(u_ode, mesh, equations, solver)
    du = wrap_array(du_ode, mesh, equations, solver)

    eta = view(u, 1, :)
    v = view(u, 2, :)
    D = view(u, 3, :)
    deta = view(du, 1, :)
    dv = view(du, 2, :)
    dD = view(du, 3, :)
    fill!(dD, zero(eltype(dD)))

    @. tmp1 = -(D * v + eta * v)
    mul!(deta, solver.D1, tmp1)

    # entropy conservative split form
    dv .= -(0.5 ./ (eta .+ D) .*
            (solver.D1 * ((eta .+ D) .* v .^ 2) + ((eta .+ D) .* v) .* (solver.D1 * v) -
             v .* (solver.D1 * ((eta .+ D) .* v))) + equations.gravity * (solver.D1 * eta))

    return nothing
end

@inline function waterheight_total(u, equations::ShallowWaterEquations1D)
    return u[1]
end

@inline function velocity(u, equations::ShallowWaterEquations1D)
    return u[2]
end

@inline function bathymetry(u, equations::ShallowWaterEquations1D)
    return -u[3]
end

@inline function waterheight(u, equations::ShallowWaterEquations1D)
    return waterheight_total(u, equations) - bathymetry(u, equations)
end

@inline function energy_total(u, equations::ShallowWaterEquations1D)
    eta, v, D = u
    e = equations.gravity * eta^2 + (D + eta) * v^2
    return e
end

@inline entropy(u, equations::ShallowWaterEquations1D) = energy_total(u, equations)

# Calculate the error for the "lake-at-rest" test case where eta should
# be a constant value over time
@inline function lake_at_rest_error(u, equations::ShallowWaterEquations1D)
    eta, _, _ = u
    return abs(equations.eta0 - eta)
end
