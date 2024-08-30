
"""
    Semidiscretization

A `struct` containing everything needed to describe a spatial semidiscretization
of an equation.
"""
struct Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                          SourceTerms, Solver, Cache}
    mesh::Mesh
    equations::Equations

    # This guy is a bit messy since we abuse it as some kind of "exact solution"
    # although this doesn't really exist...
    initial_condition::InitialCondition

    boundary_conditions::BoundaryConditions
    source_terms::SourceTerms
    solver::Solver
    cache::Cache

    function Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                                SourceTerms, Solver,
                                Cache}(mesh::Mesh, equations::Equations,
                                       initial_condition::InitialCondition,
                                       boundary_conditions::BoundaryConditions,
                                       source_terms::SourceTerms,
                                       solver::Solver,
                                       cache::Cache) where {Mesh, Equations,
                                                            InitialCondition,
                                                            BoundaryConditions, SourceTerms,
                                                            Solver, Cache}
        @assert ndims(mesh) == ndims(equations)
        @assert xmin(mesh) == xmin(solver.D1)
        @assert xmax(mesh) == xmax(solver.D1)
        @assert nnodes(mesh) == length(grid(solver))

        new(mesh, equations, initial_condition, boundary_conditions, source_terms, solver,
            cache)
    end
end

"""
    Semidiscretization(mesh, equations, initial_condition, solver;
                       source_terms=nothing,
                       boundary_conditions=boundary_condition_periodic,
                       RealT=real(solver),
                       uEltype=RealT,
                       initial_cache=(tmp1 = Array{RealT}(undef, nnodes(mesh)),))

Construct a semidiscretization of a PDE.
"""
function Semidiscretization(mesh, equations, initial_condition, solver;
                            source_terms = nothing,
                            boundary_conditions = boundary_condition_periodic,
                            # `RealT` is used as real type for node locations etc.
                            # while `uEltype` is used as element type of solutions etc.
                            RealT = real(solver), uEltype = RealT,
                            # tmp1 is needed for the `RelaxationCallback`
                            initial_cache = (tmp1 = Array{RealT}(undef, nnodes(mesh)),))
    cache = (;
             create_cache(mesh, equations, solver, initial_condition, boundary_conditions,
                          RealT, uEltype)...,
             initial_cache...)

    Semidiscretization{typeof(mesh), typeof(equations), typeof(initial_condition),
                       typeof(boundary_conditions), typeof(source_terms),
                       typeof(solver), typeof(cache)}(mesh, equations, initial_condition,
                                                      boundary_conditions, source_terms,
                                                      solver, cache)
end

function Base.show(io::IO, semi::Semidiscretization)
    @nospecialize semi # reduce precompilation time

    print(io, "Semidiscretization(")
    print(io, semi.mesh)
    print(io, ", ", semi.equations)
    print(io, ", ", semi.initial_condition)
    print(io, ", ", semi.boundary_conditions)
    print(io, ", ", semi.source_terms)
    print(io, ", ", semi.solver)
    print(io, ", cache(")
    for (idx, key) in enumerate(keys(semi.cache))
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
        println(io, "    boundary condition: ", semi.boundary_conditions)
        print(io, "    source terms: ", semi.source_terms)
    end
end

@inline Base.ndims(semi::Semidiscretization) = ndims(semi.mesh)
@inline nvariables(semi::Semidiscretization) = nvariables(semi.equations)
@inline eachvariable(semi::Semidiscretization) = eachvariable(semi.equations)
@inline nnodes(semi::Semidiscretization) = nnodes(semi.mesh)
@inline eachnode(semi::Semidiscretization) = eachnode(semi.mesh)
@inline Base.real(semi::Semidiscretization) = real(semi.solver)

"""
    grid(semi)

Get the grid of a semidiscretization.
"""
grid(semi::Semidiscretization) = grid(semi.solver)

function PolynomialBases.integrate(func, q::ArrayPartition, semi::Semidiscretization)
    integrals = zeros(real(semi), nvariables(semi))
    for v in eachvariable(semi)
        integrals[v] = integrate(func, q.x[v], semi.solver.D1)
    end
    return integrals
end
function PolynomialBases.integrate(func, quantity, semi::Semidiscretization)
    integrate(func, quantity, semi.solver.D1)
end
function PolynomialBases.integrate(q, semi::Semidiscretization)
    integrate(identity, q, semi)
end

function integrate_quantity(func, q, semi::Semidiscretization)
    quantity = zeros(eltype(q), nnodes(semi))
    integrate_quantity!(quantity, func, q, semi)
end

function integrate_quantity!(quantity, func, q, semi::Semidiscretization)
    for i in eachnode(semi)
        quantity[i] = func(get_node_vars(q, semi.equations, i), semi.equations)
    end
    integrate(quantity, semi)
end

# Obtain the function, which has an additional `!` appended to the name
inplace_version(f) = getfield(@__MODULE__, Symbol(string(nameof(f)) * "!"))

# The entropy/energy of the Svärd-Kalisch and Serre-Green-Naghdi equations
# takes the whole `q` for every point in space since it requires
# the derivative of the velocity `v_x`.
function integrate_quantity!(quantity,
                             func::Union{typeof(energy_total_modified),
                                         typeof(entropy_modified)}, q,
                             semi::Semidiscretization)
    inplace_version(func)(quantity, q, semi.equations, semi.cache)
    integrate(quantity, semi)
end

@inline function mesh_equations_solver(semi::Semidiscretization)
    @unpack mesh, equations, solver = semi
    return mesh, equations, solver
end

@inline function mesh_equations_solver_cache(semi::Semidiscretization)
    @unpack mesh, equations, solver, cache = semi
    return mesh, equations, solver, cache
end

function calc_error_norms(q, t, semi::Semidiscretization)
    calc_error_norms(q, t, semi.initial_condition, mesh_equations_solver(semi)...)
end

function rhs!(dq, q, semi::Semidiscretization, t)
    @unpack mesh, equations, initial_condition, boundary_conditions, solver, source_terms, cache = semi

    @trixi_timeit timer() "rhs!" rhs!(dq, q, t, mesh, equations, initial_condition,
                                      boundary_conditions, source_terms, solver, cache)

    return nothing
end

function compute_coefficients(func, t, semi::Semidiscretization)
    @unpack mesh, equations, solver = semi
    q = allocate_coefficients(mesh_equations_solver(semi)...)
    compute_coefficients!(q, func, t, semi)
    return q
end

function compute_coefficients!(q, func, t, semi::Semidiscretization)
    # Call `compute_coefficients` defined by the solver
    mesh, equations, solver = mesh_equations_solver(semi)
    compute_coefficients!(q, func, t, mesh,
                          is_hyperbolic_appproximation(equations), equations,
                          solver)
end

function check_bathymetry(equations, q0)
    if equations.bathymetry_type isa BathymetryFlat
        _, _, D = q0.x
        value = first(D)
        if !all(==(value), D)
            throw(ArgumentError("If the bathymetry is flat, the bathymetry should be constant."))
        end
    end
end

"""
    semidiscretize(semi::Semidiscretization, tspan)

Wrap the semidiscretization `semi` as an ODE problem in the time interval `tspan`
that can be passed to `solve` from the [SciML ecosystem](https://diffeq.sciml.ai/latest/).
"""
function semidiscretize(semi::Semidiscretization, tspan)
    q0 = compute_coefficients(semi.initial_condition, first(tspan), semi)
    check_bathymetry(semi.equations, q0)
    iip = true # is-inplace, i.e., we modify a vector when calling rhs!
    return ODEProblem{iip}(rhs!, q0, tspan, semi)
end
