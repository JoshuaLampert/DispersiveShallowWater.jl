using OrdinaryDiffEq
using SummationByPartsOperators: legendre_derivative_operator, UniformPeriodicMesh1D, couple_discontinuously
using SparseArrays: sparse
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

# create solver
p = 3 # N needs to be divisible by p + 1
Dop = legendre_derivative_operator(-1.0, 1.0, p+1)
sbp_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, N÷(p+1))
D = couple_discontinuously(Dop, sbp_mesh)
D_pl = couple_discontinuously(Dop, sbp_mesh, Val(:plus))
D_min = couple_discontinuously(Dop, sbp_mesh, Val(:minus))
D2 = sparse(D_pl) * sparse(D_min)
solver = Solver(D, D2)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation

tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)

sol = solve(ode, RK4(), save_everystep=true)

###############################################################################
# Plot numerical and analytical solution
using Plots
x = DispersiveShallowWater.grid(semi)
anim = @animate for step in 1:length(sol.u)
  t = sol.t[step]
  plot(x, sol.u[step][1, :], legend=true, ylim=(-4.0, 2.0), label="approximation", title="time t=$(round(t, digits=5))")
  analytical_sol = DispersiveShallowWater.compute_coefficients(initial_condition, solver, t, equations, mesh)
  plot!(x, analytical_sol[1, :], legend=true, label="analytical")
end
isdir("out") || mkdir("out")
gif(anim, "out/solution.gif", fps = 50)
