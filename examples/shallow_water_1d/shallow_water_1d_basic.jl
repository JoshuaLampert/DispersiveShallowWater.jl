using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: legendre_derivative_operator, UniformPeriodicMesh1D,
                                 couple_discontinuously

###############################################################################
# Semidiscretization of the shallow water equations

equations = ShallowWaterEquations1D(gravity_constant = 9.81)

initial_condition = initial_condition_radial_dambreak
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
p = 3 # N needs to be divisible by p + 1
Dop = legendre_derivative_operator(-1.0, 1.0, p + 1)
sbp_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p + 1))
D1 = couple_discontinuously(Dop, sbp_mesh)
solver = Solver(D1)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 momentum, entropy))
callbacks = CallbackSet(analysis_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
