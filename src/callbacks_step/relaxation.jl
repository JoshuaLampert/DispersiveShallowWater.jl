"""
    RelaxationCallback(invariant)

Use a relaxation method in time in order to exactly preserve the (nonlinear)
`invariant` for a conservative semidiscretization. A possible choice for
`invariant` is `invariant = entropy`.

Reference
- Hendrik Ranocha, Mohammed Sayyari, Lisandro Dalcin, Matteo Parsani, David I. Ketcheson (2019)
  Relaxation Runge–Kutta Methods: Fully-Discrete Explicit Entropy-Stable Schemes for the
  Compressible Euler and Navier–Stokes Equations
  [DOI: 10.1137/19M1263480](https://doi.org/10.1137/19M1263480)
"""
mutable struct RelaxationCallback{Invariant}
    invariant::Invariant
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:RelaxationCallback})
    @nospecialize cb # reduce precompilation time

    relaxation_callback = cb.affect!
    @unpack invariant = relaxation_callback
    print(io, "RelaxationCallback(invariant=", string(nameof(invariant)), ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:RelaxationCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        relaxation_callback = cb.affect!

        println(io, "RelaxationCallback")
        print(io, "    invariant: ", string(nameof(relaxation_callback.invariant)))
    end
end

function RelaxationCallback(; invariant)
    relaxation_callback = RelaxationCallback(invariant)

    DiscreteCallback(relaxation_callback, relaxation_callback, # the first one is the condition, the second the affect!
                     save_positions = (false, false),
                     initialize = initialize!)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: RelaxationCallback}
    return nothing
end

# this method is called to determine whether the callback should be activated
function (relaxation_callback::RelaxationCallback)(u, t, integrator)
    return true
end

# This method is called as callback during the time integration.
@inline function (relaxation_callback::RelaxationCallback)(integrator)
    semi = integrator.p
    told = integrator.tprev
    uold_ode = integrator.uprev
    uold = wrap_array(uold_ode, semi)
    tnew = integrator.t
    unew_ode = integrator.u
    unew = wrap_array(unew_ode, semi)

    gamma = one(tnew)
    terminate_integration = false
    gamma_lo = one(gamma) / 2
    gamma_hi = 3 * one(gamma) / 2

    function relaxation_functional(u, semi)
        @unpack tmp1 = semi.cache
        # modified entropy from Svärd-Kalisch equations need to take the whole vector `u` for every point in space
        if relaxation_callback.invariant isa Union{typeof(energy_total_modified), typeof(entropy_modified)}
            return integrate_quantity!(tmp1, relaxation_callback.invariant, u, semi; wrap = false)
        else
            return integrate_quantity!(tmp1, u_ode -> relaxation_callback.invariant(u_ode, semi.equations), u, semi; wrap = false)
        end
    end

    function convex_combination(gamma, uold, unew)
        @. uold + gamma * (unew - uold)
    end
    energy_old = relaxation_functional(uold, semi)

    if (relaxation_functional(convex_combination(gamma_lo, uold, unew), semi) - energy_old) *
       (relaxation_functional(convex_combination(gamma_hi, uold, unew), semi) - energy_old) >
       0
        terminate_integration = true
    else
        gamma = find_zero(g -> relaxation_functional(convex_combination(g, uold, unew),
                                                     semi) -
                               energy_old, (gamma_lo, gamma_hi), AlefeldPotraShi())
    end

    if gamma < eps(typeof(gamma))
        terminate_integration = true
    end

    unew .= convex_combination(gamma, uold, unew)
    unew_ode = reshape(unew, prod(size(unew)))
    DiffEqBase.set_u!(integrator, unew_ode)
    if !isapprox(tnew, first(integrator.opts.tstops))
        tgamma = convex_combination(gamma, told, tnew)
        DiffEqBase.set_t!(integrator, tgamma)
    end

    if terminate_integration
        terminate!(integrator)
    end

    return nothing
end
