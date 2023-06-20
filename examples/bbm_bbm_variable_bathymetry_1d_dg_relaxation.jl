using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: legendre_derivative_operator, UniformPeriodicMesh1D,
                                 couple_discontinuously
using SparseArrays: sparse

###############################################################################
# Semidiscretization of the BBM-BBM equations

equations = BBMBBMVariableEquations1D(gravity_constant = 9.81)

# initial_condition_variable_bathymetry needs periodic boundary conditions
initial_condition = initial_condition_sin_bathymetry
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver
solver = UpwindSolver(mesh, 4)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy))
relaxation_callback = RelaxationCallback(invariant = entropy)
# Always put relaxation_callback before analysis_callback to guarantee conservation of the invariant
callbacks = CallbackSet(relaxation_callback, analysis_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
