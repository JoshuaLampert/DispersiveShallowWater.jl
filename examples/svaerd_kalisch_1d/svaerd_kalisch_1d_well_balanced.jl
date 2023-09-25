using OrdinaryDiffEq
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the Svärd-Kalisch equations

equations = SvaerdKalischEquations1D(gravity_constant = 1.0, eta0 = 2.0,
                                     alpha = 0.0004040404040404049,
                                     beta = 0.49292929292929294,
                                     gamma = 0.15707070707070708)

# Setup a truly discontinuous bottom topography function for this academic
# testcase of well-balancedness. The errors from the analysis callback are
# not important but the error for this lake-at-rest test case
# `∑|eta0-eta|` should be around machine roundoff.
function initial_condition_discontinuous_well_balancedness(x, t,
                                                           equations::SvaerdKalischEquations1D,
                                                           mesh)
    # Set the background values
    eta = equations.eta0
    v = 0.0
    D = -1.0

    # Setup a discontinuous bottom topography
    if x >= 0.5 && x <= 0.75
        D = -1.5 - 0.5 * sinpi(2.0 * x)
    end

    return SVector(eta, v, D)
end

initial_condition = initial_condition_discontinuous_well_balancedness
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 200
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
solver = Solver(mesh, 4)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 momentum, entropy,
                                                                 lake_at_rest_error))
callbacks = CallbackSet(analysis_callback)

# Need a very small time step for stability
dt = 0.0002
saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), dt = dt, adaptive = false, save_everystep = false,
            callback = callbacks, saveat = saveat)
