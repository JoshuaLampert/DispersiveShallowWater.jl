using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: periodic_derivative_operator
using SparseArrays: sparse

###############################################################################
# Semidiscretization of the BBM-BBM equations

bathymetry_type = bathymetry_variable
equations = BBMBBMEquations1D(bathymetry_type, gravity_constant = 1.0, eta0 = 2.0)

initial_condition = initial_condition_discontinuous_well_balancedness
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators
accuracy_order = 4
D1 = periodic_derivative_operator(1, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
D2 = sparse(D1)^2
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
                                                                 velocity, entropy,
                                                                 lake_at_rest_error))
callbacks = CallbackSet(analysis_callback, summary_callback)

dt = 0.5
saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), dt = dt, adaptive = false, save_everystep = false,
            callback = callbacks, saveat = saveat)
