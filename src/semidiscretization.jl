
"""
    Semidiscretization

A struct containing everything needed to describe a spatial semidiscretization
of an equation.
"""
struct Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                          Solver, Cache} 

  mesh::Mesh
  equations::Equations

  # This guy is a bit messy since we abuse it as some kind of "exact solution"
  # although this doesn't really exist...
  initial_condition::InitialCondition

  boundary_conditions::BoundaryConditions
  solver::Solver
  cache::Cache

  function Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions, Solver, Cache}(
      mesh::Mesh, equations::Equations, initial_condition::InitialCondition,
      boundary_conditions::BoundaryConditions, solver::Solver,
      cache::Cache) where {Mesh, Equations, InitialCondition, BoundaryConditions, Solver, Cache}
    @assert ndims(mesh) == ndims(equations)
    @assert xmin(mesh) == xmin(solver.D)
    @assert xmax(mesh) == xmax(solver.D)
    @assert nnodes(mesh) == length(grid(solver))

    new(mesh, equations, initial_condition, boundary_conditions, solver, cache)
  end
end 


"""
    Semidiscretization(mesh, equations, initial_condition, solver;
                       boundary_conditions=boundary_condition_periodic,
                       RealT=real(solver),
                       uEltype=RealT,
                       initial_cache=NamedTuple())

Construct a semidiscretization of a PDE.
"""
function Semidiscretization(mesh, equations, initial_condition, solver;
                            boundary_conditions=boundary_condition_periodic,
                            # `RealT` is used as real type for node locations etc.
                            # while `uEltype` is used as element type of solutions etc.
                            RealT=real(solver), uEltype=RealT,
                            initial_cache=NamedTuple())

  cache = (; create_cache(mesh, equations, solver, RealT, uEltype)..., initial_cache...)

  Semidiscretization{typeof(mesh), typeof(equations), typeof(initial_condition), typeof(boundary_conditions), typeof(solver), typeof(cache)}(
    mesh, equations, initial_condition, boundary_conditions, solver, cache)
end


function Base.show(io::IO, semi::Semidiscretization)
  @nospecialize semi # reduce precompilation time

  print(io, "Semidiscretization(")
  print(io,       semi.mesh)
  print(io, ", ", semi.equations)
  print(io, ", ", semi.initial_condition)
  print(io, ", ", semi.boundary_conditions)
  print(io, ", ", semi.solver)
  print(io, ", cache(")
  for (idx,key) in enumerate(keys(semi.cache))
    idx > 1 && print(io, " ")
    print(io, key)
  end
  print(io, "))")
end

function Base.show(io::IO, ::MIME"text/plain", semi::Semidiscretization)
  @nospecialize semi # reduce precompilation time

  if get(io, :compact, false)
    show(io, semi)
  else
    println(io, "Semidiscretization")
    println(io, "    #spatial dimensions: ", ndims(semi))
    println(io, "    mesh: ", semi.mesh)
    println(io, "    equations: ", get_name(semi.equations))
    println(io, "    initial condition: ", semi.initial_condition)
    print(io, "    boundary condition: ", semi.boundary_conditions)
  end
end

@inline Base.ndims(semi::Semidiscretization) = ndims(semi.mesh)

@inline nvariables(semi::Semidiscretization) = nvariables(semi.equations)

@inline Base.real(semi::Semidiscretization) = real(semi.solver)

"""
    grid(semi)

Get the grid of a semidiscretization.
"""
grid(semi::Semidiscretization) = grid(semi.solver)

function rhs!(du_ode, u_ode, semi::Semidiscretization, t)
  @unpack mesh, equations, initial_condition, boundary_conditions, solver, cache = semi

  rhs!(du_ode, u_ode, t, mesh, equations, initial_condition, boundary_conditions, solver, cache)

  return nothing
end

@inline function set_node_vars!(u, u_node, equations, indices...)
  for v in eachvariable(equations)
    u[v, indices...] = u_node[v]
  end
  return nothing
end

function compute_coefficients(func, solver, t, equations, mesh)
  u_ode = zeros(real(solver), (nvariables(equations), nnodes(mesh)))
  compute_coefficients!(u_ode, func, solver, t, equations, mesh)
  return u_ode
end

function compute_coefficients!(u, func, solver, t, equations, mesh)
  x = grid(solver)
  for i in eachnode(solver)
    u_node = func(x[i], t, equations, mesh)
    set_node_vars!(u, u_node, equations, i)
  end
end

"""
    semidiscretize(semi::Semidiscretization, tspan)
Wrap the semidiscretization `semi` as an ODE problem in the time interval `tspan`
that can be passed to `solve` from the [SciML ecosystem](https://diffeq.sciml.ai/latest/).
"""
function semidiscretize(semi::Semidiscretization, tspan)
  u0_ode = compute_coefficients(semi.initial_condition, semi.solver, first(tspan), semi.equations, semi.mesh)
  iip = true # is-inplace, i.e., we modify a vector when calling rhs!
  return ODEProblem{iip}(rhs!, u0_ode, tspan, semi)
end
