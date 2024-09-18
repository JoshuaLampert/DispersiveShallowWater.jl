using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: fourier_derivative_operator

###############################################################################
# Semidiscretization of the BBM equation

equations = BBMEquation1D(gravity_constant = 9.81)

initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -90.0
coordinates_max = 90.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with Fourier SBP operators
D1 = fourier_derivative_operator(xmin(mesh), xmax(mesh), nnodes(mesh))
D2 = D1^2
solver = Solver(D1, D2)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified,
                                                                 hamiltonian))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
