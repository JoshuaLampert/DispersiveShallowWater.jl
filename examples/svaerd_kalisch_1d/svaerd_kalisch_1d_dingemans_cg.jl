using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: legendre_derivative_operator, UniformPeriodicMesh1D,
                                 couple_continuously, legendre_second_derivative_operator

###############################################################################
# Semidiscretization of the Svärd-Kalisch equations

equations = SvaerdKalischEquations1D(gravity_constant = 9.81, eta0 = 0.8, alpha = 0.0,
                                     beta = 0.27946992481203003, gamma = 0.0521077694235589)

initial_condition = initial_condition_dingemans
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -138.0
coordinates_max = 46.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
p = 4 # N needs to be divisible by p
D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p))
D1 = couple_continuously(D_legendre, uniform_mesh)
D2_legendre = legendre_second_derivative_operator(-1.0, 1.0, p + 1)
D2 = couple_continuously(D2_legendre, uniform_mesh)
solver = Solver(D1, D2)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 70.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy,
                                                                 entropy_modified))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 500)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
