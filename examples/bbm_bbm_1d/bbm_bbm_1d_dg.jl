using OrdinaryDiffEqTsit5
using SummationByPartsOperators: legendre_derivative_operator,
                                 UniformPeriodicMesh1D,
                                 couple_discontinuously
using SparseArrays: sparse
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the BBM-BBM equations

# or bathymetry_variable instead of bathymetry_flat
equations = BBMBBMEquations1D(bathymetry_type = bathymetry_flat, gravity_constant = 1.0)

# initial_condition_convergence_test needs periodic boundary conditions
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -35.0
coordinates_max = 35.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver
p = 3 # N needs to be divisible by p + 1
D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p + 1))
D1 = couple_discontinuously(D_legendre, uniform_mesh)
D_pl = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
D_min = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
D2 = sparse(D_pl) * sparse(D_min)
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
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
