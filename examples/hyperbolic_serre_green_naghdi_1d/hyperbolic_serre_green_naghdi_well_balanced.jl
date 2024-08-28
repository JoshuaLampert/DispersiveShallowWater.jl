using OrdinaryDiffEq
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope,
                                                  lambda = 500.0,
                                                  gravity_constant = 1.0,
                                                  eta0 = 2.0)

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
