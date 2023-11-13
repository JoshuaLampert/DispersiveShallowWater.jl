"""
    AbstractSolver

An abstract supertype of specific solvers.
"""
abstract type AbstractSolver end

"""
    Solver

A `struct` that holds the summation by parts (SBP) operators that are used for the spatial discretization.
"""
struct Solver{RealT <: Real} <: AbstractSolver
    D1::AbstractDerivativeOperator{RealT}
    D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}}

    function Solver{RealT}(D1::AbstractDerivativeOperator{RealT},
                           D2::Union{AbstractDerivativeOperator{RealT},
                                     AbstractMatrix{RealT}}) where {RealT}
        @assert derivative_order(D1) == 1
        if D2 isa AbstractDerivativeOperator
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

Create a solver, where `D1` is an `AbstractDerivativeOperator`  from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
of first `derivative_order` and `D2` is an `AbstractDerivativeOperator` of second `derivative_order` or an `AbstractMatrix`.
Both summation by parts operators should be associated with the same grid.
"""
function Solver(D1::AbstractDerivativeOperator{RealT},
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
    q = wrap_array(u_ode, mesh, equations, solver)
    q_exact = zeros(real(solver), (nvariables(equations), nnodes(mesh)))
    for i in eachnode(solver)
        q_exact[:, i] = initial_condition(x[i], t, equations, mesh)
    end
    l2_error = zeros(real(solver), nvariables(equations))
    linf_error = similar(l2_error)
    for v in eachvariable(equations)
        @views diff = q[v, :] - q_exact[v, :]
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
        local_source = source_terms(view(q, :, i), x[i], t, equations)
        for v in eachvariable(equations)
            dq[v, i] += local_source[v]
        end
    end
    return nothing
end
