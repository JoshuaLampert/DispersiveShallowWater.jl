using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: fourier_derivative_operator

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

bathymetry_type = bathymetry_flat
equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type;
                                                  lambda = 500.0,
                                                  gravity_constant = 9.81)

initial_condition = initial_condition_soliton
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with Fourier pseudospectral collocation method
D1 = fourier_derivative_operator(xmin(mesh), xmax(mesh), nnodes(mesh))
solver = Solver(D1, nothing)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * gravity_constant(equations))) # one period
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified))
relaxation_callback = RelaxationCallback(invariant = entropy_modified)
# Always put relaxation_callback before analysis_callback to guarantee conservation of the invariant
callbacks = CallbackSet(relaxation_callback, analysis_callback, summary_callback)

# optimized time integration methods like this one are much more efficient
# for stiff problems (λ big) than standard methods like Tsit5()
alg = RDPK3SpFSAL35()
sol = solve(ode, alg, abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks)
