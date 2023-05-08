using OrdinaryDiffEq
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the BBM-BBM equations

equations = BBMBBMEquations1D(gravity_constant=1.0, D=1.0)

# initial_condition_convergence_test needs periodic boundary conditions
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -35.0
coordinates_max = 35.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
solver = Solver(mesh, 4)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation

tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)

sol = solve(ode, RK4(), save_everystep=true)

