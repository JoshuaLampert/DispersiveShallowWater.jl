using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
using SparseArrays: sparse

###############################################################################
# Semidiscretization of the BBM-BBM equations

# or bathymetry_variable instead of bathymetry_flat
equations = BBMBBMEquations1D(bathymetry_type = bathymetry_flat, gravity = 9.81)

# initial_condition_convergence_test needs periodic boundary conditions
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -35.0
coordinates_max = 35.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver
accuracy_order = 4
D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                      accuracy_order = accuracy_order, xmin = mesh.xmin, xmax = mesh.xmax,
                      N = mesh.N)
D2 = sparse(D1.plus) * sparse(D1.minus)
solver = Solver(D1, D2)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy))
relaxation_callback = RelaxationCallback(invariant = entropy)
# Always put relaxation_callback before analysis_callback to guarantee conservation of the invariant
callbacks = CallbackSet(relaxation_callback, analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
