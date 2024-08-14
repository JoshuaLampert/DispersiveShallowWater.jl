"""
    AbstractSolver

An abstract supertype of specific solvers.
"""
abstract type AbstractSolver end

"""
    Solver

A `struct` that holds the summation by parts (SBP) operators that are used for the spatial discretization.
"""
struct Solver{RealT <: Real,
              FirstDerivative <: AbstractDerivativeOperator{RealT},
              SecondDerivative <:
              Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}, Nothing}} <:
       AbstractSolver
    D1::FirstDerivative
    D2::SecondDerivative

    function Solver{RealT, FirstDerivative, SecondDerivative}(D1::FirstDerivative,
                                                              D2::SecondDerivative) where {
                                                                                           RealT,
                                                                                           FirstDerivative,
                                                                                           SecondDerivative
                                                                                           }
        @assert derivative_order(D1) == 1
        if D2 isa AbstractDerivativeOperator &&
           !(D2 isa SummationByPartsOperators.FourierPolynomialDerivativeOperator)
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
    Solver{real(D1), typeof(D1), typeof(D2)}(D1, D2)
end

# Also allow to pass custom SBP operators (for convenience without explicitly specifying the type)
"""
    Solver(D1, D2)

Create a solver, where `D1` is an `AbstractDerivativeOperator`
from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
of first `derivative_order` and `D2` is an `AbstractDerivativeOperator`
of second `derivative_order` or an `AbstractMatrix`. It can also be `nothing`
if no second derivative is used by the discretization.
Both summation-by-parts operators should be associated with the same grid.
"""
function Solver(D1::AbstractDerivativeOperator{RealT},
                D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT},
                          Nothing}) where {RealT}
    Solver{RealT, typeof(D1), typeof(D2)}(D1, D2)
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
@inline Base.real(::Solver{RealT}) where {RealT} = RealT

# Adapted from Trixi.jl
# https://github.com/trixi-framework/Trixi.jl/blob/75d8c67629562efd24b2a04e46d22b0a1f4f572c/src/solvers/dg.jl#L539
@inline function get_node_vars(q, equations, indices...)
    # There is a cut-off at `n == 10` inside of the method
    # `ntuple(f::F, n::Integer) where F` in Base at ntuple.jl:17
    # in Julia `v1.5`, leading to type instabilities if
    # more than ten variables are used. That's why we use
    # `Val(...)` below.
    # We use `@inline` to make sure that the `getindex` calls are
    # really inlined, which might be the default choice of the Julia
    # compiler for standard `Array`s but not necessarily for more
    # advanced array types such as `PtrArray`s, cf.
    # https://github.com/JuliaSIMD/VectorizationBase.jl/issues/55
    SVector(ntuple(@inline(v->q.x[v][indices...]), Val(nvariables(equations))))
end

@inline function set_node_vars!(q, q_node, equations, indices...)
    for v in eachvariable(equations)
        q.x[v][indices...] = q_node[v]
    end
    return nothing
end

function allocate_coefficients(mesh::Mesh1D, equations, solver::AbstractSolver)
    return ArrayPartition(ntuple(_ -> zeros(real(solver), nnodes(mesh)),
                                 Val(nvariables(equations))))
end

function compute_coefficients!(q, func, t, mesh::Mesh1D, equations, solver::AbstractSolver)
    x = grid(solver)
    for i in eachnode(solver)
        q_node = func(x[i], t, equations, mesh)
        set_node_vars!(q, q_node, equations, i)
    end
end

function calc_error_norms(q, t, initial_condition, mesh::Mesh1D, equations,
                          solver::AbstractSolver)
    x = grid(solver)
    q_exact = zeros(real(solver), (nvariables(equations), nnodes(mesh)))
    for i in eachnode(solver)
        q_exact[:, i] = initial_condition(x[i], t, equations, mesh)
    end
    l2_error = zeros(real(solver), nvariables(equations))
    linf_error = similar(l2_error)
    for v in eachvariable(equations)
        @views diff = q.x[v] - q_exact[v, :]
        l2_error[v] = integrate(q -> q^2, diff, solver.D1) |> sqrt
        linf_error[v] = maximum(abs.(diff))
    end
    return l2_error, linf_error
end

function calc_sources!(dq, q, t, source_terms::Nothing,
                       equations::AbstractEquations{1}, solver::Solver)
    return nothing
end

function calc_sources!(dq, q, t, source_terms,
                       equations::AbstractEquations{1}, solver::Solver)
    x = grid(solver)
    for i in eachnode(solver)
        local_source = source_terms(get_node_vars(q, equations, i), x[i], t, equations)
        for v in eachvariable(equations)
            dq.x[v][i] += local_source[v]
        end
    end
    return nothing
end
