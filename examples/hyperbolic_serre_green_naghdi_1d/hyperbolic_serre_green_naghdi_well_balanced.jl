using OrdinaryDiffEq
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

bathymetry_type = bathymetry_mild_slope
equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type;
                                                  lambda = 500.0,
                                                  gravity_constant = 1.0,
                                                  eta0 = 2.0)

# Setup a truly discontinuous bottom topography function for this academic
# testcase of well-balancedness. The errors from the analysis callback are
# not important but the error for this lake-at-rest test case
# `∫|η-η₀|` should be around machine roundoff.
function initial_condition_discontinuous_well_balancedness(x, t,
                                                           equations::HyperbolicSerreGreenNaghdiEquations1D,
                                                           mesh)
    # Set the background values
    eta = equations.eta0
    v = 0.0
    D = equations.eta0 - 1.0

    # Setup a discontinuous bottom topography
    if x >= 0.5 && x <= 0.75
        D = equations.eta0 - 1.5 - 0.5 * sinpi(2.0 * x)
    end

    w = 0.0 # -h v_x
    H = eta + D - equations.eta0

    return SVector(eta, v, D, w, H)
end

initial_condition = initial_condition_discontinuous_well_balancedness
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
accuracy_order = 4
solver = Solver(mesh, accuracy_order)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 momentum,
                                                                 entropy_modified,
                                                                 lake_at_rest_error))
callbacks = CallbackSet(analysis_callback, summary_callback)

# optimized time integration methods like this one are much more efficient
# for stiff problems (λ big) than standard methods like Tsit5()
alg = RDPK3SpFSAL35()
sol = solve(ode, alg, dt = 2.5e-4, adaptive = false, callback = callbacks)
