"""
    AnalysisCallback(semi; interval=0,
                           extra_analysis_errors=Symbol[],
                           extra_analysis_integrals=())

Analyze a numerical solution every `interval` time steps.
The L2- and the L∞-norm for each component are computed by default.
Additional errors can be computed, e.g. by passing `extra_analysis_errors = (:conservation_error,)`.

Further scalar functions `func` in `extra_analysis_integrals` are applied to the numerical
solution and integrated over the computational domain. Some examples for this are
[`entropy`](@ref), and [`energy_total`](@ref).
You can also write your own function with the same signature as the examples listed above and
pass it via `extra_analysis_integrals`.
The computed errors and intergrals are saved for each timestep and can be obtained by calling
[`errors`](@ref) and [`integrals`](@ref).
"""
mutable struct AnalysisCallback{T, AnalysisIntegrals, InitialStateIntegrals}
  interval::Int
  analysis_errors::Vector{Symbol}
  analysis_integrals::AnalysisIntegrals
  initial_state_integrals::InitialStateIntegrals
  tstops::Vector{Float64}
  errors::Vector{Matrix{T}}
  integrals::Vector{Vector{T}}
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:AnalysisCallback})
  @nospecialize cb # reduce precompilation time

  analysis_callback = cb.affect!
  @unpack interval = analysis_callback
  print(io, "AnalysisCallback(interval=", interval, ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:AnalysisCallback})
  @nospecialize cb # reduce precompilation time

  if get(io, :compact, false)
    show(io, cb)
  else
    analysis_callback = cb.affect!
    @unpack interval, analysis_errors, analysis_integrals = analysis_callback

    println(io, "AnalysisCallback")
    println(io, "    interval: ", interval)
    for (idx, error) in enumerate(analysis_errors)
      println(io, "    error ", idx, ": ", error)
    end
    for (idx, integral) in enumerate(analysis_integrals)
      println(io, "    integral ", idx, ": ", integral)
    end
  end
end

function AnalysisCallback(semi::Semidiscretization; kwargs...)
  mesh, equations, solver = mesh_equations_solver(semi)
  AnalysisCallback(mesh, equations, solver; kwargs...)
end

function AnalysisCallback(mesh, equations::AbstractEquations, solver;
                          interval = 0,
                          extra_analysis_errors = Symbol[],
                          analysis_errors = union(default_analysis_errors(equations),
                                                  extra_analysis_errors),
                          extra_analysis_integrals = (),
                          analysis_integrals = union(default_analysis_integrals(equations),
                                                     extra_analysis_integrals))
  # Decide when the callback is activated.
  # With error-based step size control, some steps can be rejected. Thus,
  #   `integrator.iter >= integrator.stats.naccept`
  #    (total #steps)       (#accepted steps)
  # We need to check the number of accepted steps since callbacks are not
  # activated after a rejected step.
  condition = (u, t, integrator) -> interval > 0 &&
    ((integrator.stats.naccept % interval == 0 &&
      !(integrator.stats.naccept == 0 && integrator.iter > 0)) ||
     isfinished(integrator))
  for extra_analysis_error in extra_analysis_errors
    if extra_analysis_error != :conservation_error
      @warn "Extra analysis error $extra_analysis_error is not supported and will be ignored."
    end
  end
  analysis_callback = AnalysisCallback(interval, analysis_errors, Tuple(analysis_integrals),
                                       SVector(ntuple(_ -> zero(real(solver)),
                                                      Val(nvariables(equations)))),
                                       Vector{Float64}(),
                                       Vector{Matrix{real(solver)}}(),
                                       Vector{Vector{real(solver)}}())

  DiscreteCallback(condition, analysis_callback,
                   save_positions = (false, false),
                   initialize = initialize!)
end

"""
    tstops(analysis_callback)

Return the time values that correspond to the saved values of the [`errors`](@ref) and [`integrals`](@ref).
"""
function tstops(cb::DiscreteCallback{
                                     Condition,
                                     Affect!
                                     }) where {Condition,
                                               Affect! <:
                                               AnalysisCallback}
  analysis_callback = cb.affect!
  return analysis_callback.tstops
end

"""
    errors(analysis_callback)
    
Return the computed errors for each timestep. The shape is (nerrors, nvariables, ntimesteps).
"""
function errors(cb::DiscreteCallback{
                                     Condition,
                                     Affect!
                                     }) where {Condition,
                                               Affect! <:
                                               AnalysisCallback}
  analysis_callback = cb.affect!
  # reshape to return 3D array of shape (nerrors, nvariables, ntimesteps)
  return reshape(reduce(hcat, analysis_callback.errors),
                 size(analysis_callback.errors[1])..., :)
end

"""
    integrals(analysis_callback)

Return the computed integrals for each timestep. The shape is (nintegrals, ntimesteps).
"""
function integrals(cb::DiscreteCallback{
                                        Condition,
                                        Affect!
                                        }) where {Condition,
                                                  Affect! <:
                                                  AnalysisCallback}
  analysis_callback = cb.affect!
  # reshape to return Matrix of shape (nintegrals, ntimesteps)
  return reduce(hcat, analysis_callback.integrals)
end

"""
    integral_names(analysis_callback)

Return the names of the computed integrals.
"""
function integral_names(cb::DiscreteCallback{
                                             Condition,
                                             Affect!
                                             }) where {Condition,
                                                       Affect! <:
                                                       AnalysisCallback}
  analysis_callback = cb.affect!
  return nameof.(analysis_callback.analysis_integrals)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, t,
                     integrator) where {Condition, Affect! <: AnalysisCallback}
  semi = integrator.p
  initial_state_integrals = integrate(u_ode, semi)

  analysis_callback = cb.affect!
  analysis_callback.initial_state_integrals = initial_state_integrals

  analysis_callback(integrator)
  return nothing
end

function (analysis_callback::AnalysisCallback)(integrator)
  semi = integrator.p
  mesh, equations, solver = mesh_equations_solver(semi)
  @unpack t = integrator

  # Calculate current time derivative (needed for semidiscrete entropy time derivative, residual, etc.)
  du_ode = first(get_tmp_cache(integrator))
  # `integrator.f` is usually just a call to `rhs!`
  # However, we want to allow users to modify the ODE RHS outside of Trixi.jl
  # and allow us to pass a combined ODE RHS to OrdinaryDiffEq, e.g., for
  # hyperbolic-parabolic systems.
  integrator.f(du_ode, integrator.u, semi, t)
  l2_error, linf_error = analysis_callback(integrator.u, t, semi)

  # avoid re-evaluating possible FSAL stages
  u_modified!(integrator, false)

  # Return errors for EOC analysis
  return l2_error, linf_error
end

# This method is just called internally from `(analysis_callback::AnalysisCallback)(integrator)`
# and serves as a function barrier. Additionally, it makes the code easier to profile and optimize.
function (analysis_callback::AnalysisCallback)(u_ode, t, semi)
  _, equations, solver = mesh_equations_solver(semi)
  @unpack analysis_errors, analysis_integrals, tstops, errors, integrals = analysis_callback

  push!(tstops, t)
  # Calculate L2/Linf errors, which are also returned
  l2_error, linf_error = calc_error_norms(u_ode, t, semi)
  current_errors = zeros(real(semi), (length(analysis_errors), nvariables(equations)))
  current_errors[1, :] = l2_error
  current_errors[2, :] = linf_error

  # Conservation error
  if :conservation_error in analysis_errors
    @unpack initial_state_integrals = analysis_callback
    state_integrals = integrate(u_ode, semi)
    current_errors[3, :] = abs.(state_integrals - initial_state_integrals)
  end
  push!(errors, current_errors)

  # additional integrals
  current_integrals = zeros(real(semi), length(analysis_integrals))
  analyze_integrals!(current_integrals, 1, analysis_integrals, u_ode, t, semi)
  push!(integrals, current_integrals)

  return l2_error, linf_error
end

# Iterate over tuples of analysis integrals in a type-stable way using "lispy tuple programming".
function analyze_integrals!(current_integrals, i, analysis_integrals::NTuple{N, Any}, u_ode,
                            t, semi) where {N}

  # Extract the first analysis integral and process it; keep the remaining to be processed later
  quantity = first(analysis_integrals)
  remaining_quantities = Base.tail(analysis_integrals)

  res = analyze(quantity, u_ode, t, semi)
  current_integrals[i] = res

  # Recursively call this method with the unprocessed integrals
  analyze_integrals!(current_integrals, i + 1, remaining_quantities, u_ode, t, semi)
  return nothing
end

# terminate the type-stable iteration over tuples
function analyze_integrals!(current_integrals, i, analysis_integrals::Tuple{}, u_ode, t,
                            semi)
  nothing
end

function analyze(quantity, u_ode, t, semi::Semidiscretization)
  integrate_quantity(u_ode -> quantity(u_ode, semi.equations), u_ode, semi)
end
