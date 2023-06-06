@inline function isfinished(integrator)
  # Checking for floating point equality is OK here as `DifferentialEquations.jl`
  # sets the time exactly to the final time in the last iteration
  return integrator.t == last(integrator.sol.prob.tspan) ||
         isempty(integrator.opts.tstops) ||
         integrator.iter == integrator.opts.maxiters
end

include("relaxation.jl")
include("analysis.jl")
