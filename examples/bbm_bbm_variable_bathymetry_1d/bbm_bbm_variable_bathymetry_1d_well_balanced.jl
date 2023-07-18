using OrdinaryDiffEq
using DispersiveShallowWater
using SummationByPartsOperators: periodic_derivative_operator
using SparseArrays: sparse

###############################################################################
# Semidiscretization of the BBM-BBM equations

equations = BBMBBMVariableEquations1D(gravity_constant = 9.81, eta0 = 2.0)

# Setup a truly discontinuous bottom topography function for this academic
# testcase of well-balancedness. The errors from the analysis callback are
# not important but the error for this lake-at-rest test case
# `∑|eta0-eta|` should be around machine roundoff.
function initial_condition_discontinuous_well_balancedness(x, t,
                                                           equations::BBMBBMVariableEquations1D,
                                                           mesh)
    # Set the background values
    eta = equations.eta0
    v = 0.0
    D = -1.0

    # Setup a discontinuous bottom topography
    if x >= 0.5 && x <= 0.75
        D = -1.5 - 0.5 * sinpi(2.0 * x)
    end

    return SVector(eta, v, D)
end

initial_condition = initial_condition_discontinuous_well_balancedness
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N + 1)

# create solver with periodic SBP operators
D1 = periodic_derivative_operator(1, 4, mesh.xmin, mesh.xmax, mesh.N)
D2 = sparse(D1)^2
solver = Solver(D1, D2)

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
                                                                 velocity, entropy,
                                                                 lake_at_rest_error))
callbacks = CallbackSet(analysis_callback)

dt = 0.5
saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), dt = dt, adaptive = false, save_everystep = false,
            callback = callbacks, saveat = saveat)
