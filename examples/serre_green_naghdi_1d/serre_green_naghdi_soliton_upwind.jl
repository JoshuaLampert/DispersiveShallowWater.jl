using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
using SparseArrays: sparse

###############################################################################
# Semidiscretization of the Serre-Green-Naghdi equations

equations = SerreGreenNaghdiEquations1D(gravity_constant = 9.81)

# initial_condition_convergence_test needs periodic boundary conditions
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic upwind SBP operators of accuracy order 4
# (note that only central-type with even order are used)
accuracy_order = 4
D1 = upwind_operators(periodic_derivative_operator;
                      derivative_order = 1,
                      accuracy_order = accuracy_order - 1,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1, nothing)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * equations.gravity)) # one period
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
