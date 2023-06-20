"""
    AbstractSolver

An abstract supertype of specific solvers.
"""
abstract type AbstractSolver end

"""
    Solver

A struct that holds the summation by parts (SBP) operators that are used for the spatial discretization.
"""
struct Solver{RealT <: Real} <: AbstractSolver
  D1::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}
  D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}

  function Solver{RealT}(D1::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}
                                   },
                         D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}
                                   }) where {RealT}
    if D1 isa AbstractDerivativeOperator{RealT} && D2 isa AbstractDerivativeOperator{RealT}
      @assert grid(D1) == grid(D2)
    end
    if D1 isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D1) == 1
    end
    if D2 isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D2) == 2
    end
    new(D1, D2)
  end
end

"""
    Solver(mesh, accuracy_order)

Create a solver, where the summation by parts (SBP) operators are of order `accuracy_order` and
associated to the `mesh`.
"""
function Solver(mesh, accuracy_order)
  D1 = periodic_derivative_operator(1, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
  D2 = periodic_derivative_operator(2, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
  @assert real(D1) == real(D2)
  Solver{real(D1)}(D1, D2)
end

# Also allow to pass custom SBP operators (for convenience without explicitly specifying the type)
"""
    Solver(D1, D2)

Create a solver, where `D1` is an `AbstractDerivativeOperator` of first `derivative_order` and `D2`
is an `AbstractDerivativeOperator` of second `derivative_order` from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).
Both summation by parts operators should be associated with the same grid.
"""
function Solver(D1::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}},
                D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}) where {
                                                                                            RealT
                                                                                            }
  Solver{RealT}(D1, D2)
end

function Base.show(io::IO, solver::Solver{RealT}) where {RealT}
  print(io, "Solver{", RealT, "}")
end

function Base.show(io::IO, ::MIME"text/plain", solver::Solver{RealT}) where {RealT}
  if get(io, :compact, false)
    show(io, solver)
  else
    println(io, "Solver{", RealT, "}")
    println(io, "    D1: ", solver.D1)
    print(io, "    D2: ", solver.D2)
  end
end

grid(solver::Solver) = grid(solver.D1)

@inline eachnode(solver::Solver) = Base.OneTo(length(grid(solver)))
@inline Base.real(solver::Solver{RealT}) where {RealT} = RealT

"""
    UpwindSolver

A struct that holds the upwind summation by parts (SBP) operators that are used for the spatial discretization.
"""
struct UpwindSolver{RealT <: Real} <: AbstractSolver
  D1::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}
  D_pl::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}
  D_min::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}
  D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}

  function UpwindSolver{RealT}(D1::Union{AbstractDerivativeOperator{RealT},
                                         AbstractMatrix{RealT}},
                               D_pl::Union{AbstractDerivativeOperator{RealT},
                                           AbstractMatrix{RealT}},
                               D_min::Union{AbstractDerivativeOperator{RealT},
                                            AbstractMatrix{RealT}
                                            }) where {RealT}
    if D1 isa AbstractDerivativeOperator{RealT} &&
       D_pl isa AbstractDerivativeOperator{RealT} &&
       D_min isa AbstractDerivativeOperator{RealT}
      @assert grid(D1) == grid(D_pl) == grid(D_min)
    end
    if D1 isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D1) == 1
    end
    if D_pl isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D_pl) == 1
    end
    if D_min isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D_min) == 1
    end
    new(D1, D_pl, D_min, sparse(D_pl) * sparse(D_min))
  end
end

"""
    UpwindSolver(mesh, accuracy_order)

Create a solver, where the summation by parts (SBP) operators are of order `accuracy_order` and
associated to the `mesh`.
"""
function UpwindSolver(mesh, accuracy_order)
  Dop = legendre_derivative_operator(mesh.xmin, mesh.xmax, accuracy_order)
  sbp_mesh = UniformPeriodicMesh1D(-1.0, 1.0, div(mesh.N, accuracy_order))
  D1 = couple_discontinuously(Dop, sbp_mesh)
  D_pl = couple_discontinuously(Dop, sbp_mesh, Val(:plus))
  D_min = couple_discontinuously(Dop, sbp_mesh, Val(:minus))
  @assert real(D1) == real(D_pl) == real(D_min)
  UpwindSolver{real(D1)}(D1, D_pl, D_min)
end

# Also allow to pass custom SBP operators (for convenience without explicitly specifying the type)
"""
    UpwindSolver(D1, D_pl, D_min)

Create a solver, where `D1` is an `AbstractDerivativeOperator` of first `derivative_order`, `D_pl` is a
positive upwind `AbstractDerivativeOperator` of first `derivative_order` and `D_min` is a negative upwind
`AbstractDerivativeOperator` of first `derivative_order` from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).
All three summation by parts operators should be associated with the same grid.
"""
function UpwindSolver(D1::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}},
                      D_pl::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}},
                      D_min::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}
                                   }) where {
                                             RealT
                                             }
  UpwindSolver{RealT}(D1, D_pl, D_mi)
end

function Base.show(io::IO, solver::UpwindSolver{RealT}) where {RealT}
  print(io, "UpwindSolver{", RealT, "}")
end

function Base.show(io::IO, ::MIME"text/plain", solver::UpwindSolver{RealT}) where {RealT}
  if get(io, :compact, false)
    show(io, solver)
  else
    println(io, "UpwindSolver{", RealT, "}")
    println(io, "    D_pl: ", solver.D_pl)
    print(io, "    D_min: ", solver.D_min)
  end
end

grid(solver::UpwindSolver) = grid(solver.D_pl)

@inline eachnode(solver::UpwindSolver) = Base.OneTo(length(grid(solver)))
@inline Base.real(solver::UpwindSolver{RealT}) where {RealT} = RealT

@inline function set_node_vars!(u, u_node, equations, indices...)
  for v in eachvariable(equations)
    u[v, indices...] = u_node[v]
  end
  return nothing
end

function allocate_coefficients(mesh::Mesh1D, equations, solver::AbstractSolver)
  # cf. wrap_array
  zeros(real(solver), nvariables(equations) * nnodes(mesh)^ndims(mesh))
end

function wrap_array(u_ode, mesh, equations, solver)
  unsafe_wrap(Array{eltype(u_ode), ndims(mesh) + 1}, pointer(u_ode),
              (nvariables(equations), ntuple(_ -> nnodes(mesh), ndims(mesh))...))
end

function compute_coefficients!(u, func, t, mesh::Mesh1D, equations, solver::AbstractSolver)
  x = grid(solver)
  for i in eachnode(solver)
    u_node = func(x[i], t, equations, mesh)
    set_node_vars!(u, u_node, equations, i)
  end
end

function calc_error_norms(u_ode, t, initial_condition, mesh::Mesh1D, equations,
                          solver::AbstractSolver)
  x = grid(solver)
  u = wrap_array(u_ode, mesh, equations, solver)
  u_exact = zeros(real(solver), (nvariables(equations), nnodes(mesh)))
  for i in eachnode(solver)
    u_exact[:, i] = initial_condition(x[i], t, equations, mesh)
  end
  l2_error = zeros(real(solver), nvariables(equations))
  linf_error = similar(l2_error)
  for v in eachvariable(equations)
    @views diff = u[v, :] - u_exact[v, :]
    l2_error[v] = integrate(u -> u^2, diff, solver.D1) |> sqrt
    linf_error[v] = maximum(abs.(diff))
  end
  return l2_error, linf_error
end
