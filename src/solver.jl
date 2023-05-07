"""
    Solver

A struct that holds the summation by parts (SBP) operators that are used for the spatial discretization.
"""
struct Solver{RealT <: Real}
  D::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}
  D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}

  function Solver{RealT}(D::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}},
                         D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}) where {RealT}
    if D isa AbstractDerivativeOperator{RealT} && D2 isa AbstractDerivativeOperator{RealT}
      @assert grid(D) == grid(D2)
    end
    if D isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D) == 1
    end
    if D2 isa AbstractDerivativeOperator{RealT}
      @assert derivative_order(D2) == 2
    end
    new(D, D2)
  end
end

"""
    Solver(mesh, accuracy_order)

Create a solver, where the summation by parts (SBP) operators are of order `accuracy_order` and
associated to the `mesh`.
"""
function Solver(mesh, accuracy_order)
  D = periodic_derivative_operator(1, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
  D2 = periodic_derivative_operator(2, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
  @assert real(D) == real(D2)
  Solver{real(D)}(D, D2)
end

# Also allow to pass custom SBP operators (for convenience without explicitly specifying the type)
"""
    Solver(D, D2)

Create a solver, where `D` is an `AbstractDerivativeOperator` of first `derivative_order` and `D2`
is an `AbstractDerivativeOperator` of second `derivative_order` from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).
Both summation by parts operators should be associated with the same grid.
"""
function Solver(D::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}},
                D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}) where {RealT}
  Solver{RealT}(D, D2)
end

function Base.show(io::IO, solver::Solver{RealT}) where {RealT}
  print(io, "Solver{", RealT, "}")
end

function Base.show(io::IO, ::MIME"text/plain", solver::Solver{RealT}) where {RealT}
  if get(io, :compact, false)
    show(io, solver)
  else
    println(io, "Solver{", RealT, "}")
    println(io, "    D: ", solver.D)
    print(io, "    D2: ", solver.D2)
  end
end

grid(solver::Solver) = grid(solver.D)

@inline eachnode(solver::Solver) = Base.OneTo(length(grid(solver)))
@inline Base.real(solver::Solver{RealT}) where {RealT} = RealT
