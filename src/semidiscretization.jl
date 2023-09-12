
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

    function Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                                Solver,
                                Cache}(mesh::Mesh, equations::Equations,
                                       initial_condition::InitialCondition,
                                       boundary_conditions::BoundaryConditions,
                                       solver::Solver,
                                       cache::Cache) where {Mesh, Equations,
                                                            InitialCondition,
                                                            BoundaryConditions, Solver,
                                                            Cache}
        @assert ndims(mesh) == ndims(equations)
        @assert xmin(mesh) == xmin(solver.D1)
        @assert xmax(mesh) == xmax(solver.D1)
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
                            boundary_conditions = boundary_condition_periodic,
                            # `RealT` is used as real type for node locations etc.
                            # while `uEltype` is used as element type of solutions etc.
                            RealT = real(solver), uEltype = RealT,
                            initial_cache = NamedTuple())
    cache = (;
             create_cache(mesh, equations, solver, initial_condition, RealT, uEltype)...,
             initial_cache...)

    Semidiscretization{typeof(mesh), typeof(equations), typeof(initial_condition),
                       typeof(boundary_conditions), typeof(solver), typeof(cache)}(mesh,
                                                                                   equations,
                                                                                   initial_condition,
                                                                                   boundary_conditions,
                                                                                   solver,
                                                                                   cache)
end

function Base.show(io::IO, semi::Semidiscretization)
    @nospecialize semi # reduce precompilation time

    print(io, "Semidiscretization(")
    print(io, semi.mesh)
    print(io, ", ", semi.equations)
    print(io, ", ", semi.initial_condition)
    print(io, ", ", semi.boundary_conditions)
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
        print(io, "    boundary condition: ", semi.boundary_conditions)
    end
end

@inline Base.ndims(semi::Semidiscretization) = ndims(semi.mesh)
@inline nvariables(semi::Semidiscretization) = nvariables(semi.equations)
@inline nnodes(semi::Semidiscretization) = nnodes(semi.mesh)
@inline eachnode(semi::Semidiscretization) = eachnode(semi.mesh)
@inline Base.real(semi::Semidiscretization) = real(semi.solver)

"""
    grid(semi)

Get the grid of a semidiscretization.
"""
grid(semi::Semidiscretization) = grid(semi.solver)

function PolynomialBases.integrate(func, u_ode, semi::Semidiscretization; wrap = true)
    if wrap == true
        u = wrap_array(u_ode, semi)
        integrals = zeros(real(semi), nvariables(semi))
        for v in eachvariable(semi.equations)
            integrals[v] = integrate(func, u[v, :], semi.solver.D1)
        end
        return integrals
    else
        integrate(func, u_ode, semi.solver.D1)
    end
end
function PolynomialBases.integrate(u, semi::Semidiscretization; wrap = true)
    integrate(identity, u, semi; wrap = wrap)
end

function integrate_quantity(func, u_ode, semi::Semidiscretization; wrap = true)
    if wrap == true
        u = wrap_array(u_ode, semi)
    else
        u = u_ode
    end
    quantity = zeros(eltype(u), size(u, 2))
    for i in 1:size(u, 2)
        quantity[i] = func(view(u, :, i))
    end
    integrate(quantity, semi; wrap = false)
end

function integrate_quantity!(quantity, func, u_ode, semi::Semidiscretization; wrap = true)
    if wrap == true
        u = wrap_array(u_ode, semi)
    else
        u = u_ode
    end
    for i in 1:size(u, 2)
        quantity[i] = func(view(u, :, i))
    end
    integrate(quantity, semi; wrap = false)
end

# modified entropy from Svärd-Kalisch equations need to take the whole vector `u` for every point in space
function integrate_quantity(func::Union{typeof(energy_total_modified),
                                        typeof(entropy_modified)}, u_ode,
                            semi::Semidiscretization; wrap = true)
    if wrap == true
        u = wrap_array(u_ode, semi)
    else
        u = u_ode
    end
    quantity = func(u, semi.equations, semi.cache)
    integrate(quantity, semi; wrap = false)
end

function integrate_quantity!(quantity,
                             func::Union{typeof(energy_total_modified),
                                         typeof(entropy_modified)}, u_ode,
                             semi::Semidiscretization; wrap = true)
    if wrap == true
        u = wrap_array(u_ode, semi)
    else
        u = u_ode
    end
    quantity = func(u, semi.equations, semi.cache)
    integrate(quantity, semi; wrap = false)
end

@inline function mesh_equations_solver(semi::Semidiscretization)
    @unpack mesh, equations, solver = semi
    return mesh, equations, solver
end

@inline function mesh_equations_solver_cache(semi::Semidiscretization)
    @unpack mesh, equations, solver, cache = semi
    return mesh, equations, solver, cache
end

function calc_error_norms(u, t, semi::Semidiscretization)
    calc_error_norms(u, t, semi.initial_condition, mesh_equations_solver(semi)...)
end

function wrap_array(u_ode, semi::Semidiscretization)
    wrap_array(u_ode, mesh_equations_solver(semi)...)
end

function rhs!(du_ode, u_ode, semi::Semidiscretization, t)
    @unpack mesh, equations, initial_condition, boundary_conditions, solver, cache = semi

    rhs!(du_ode, u_ode, t, mesh, equations, initial_condition, boundary_conditions, solver,
         cache)

    return nothing
end

function compute_coefficients(func, t, semi::Semidiscretization)
    @unpack mesh, equations, solver = semi
    u_ode = allocate_coefficients(mesh_equations_solver(semi)...)
    compute_coefficients!(u_ode, func, t, semi)
    return u_ode
end

function compute_coefficients!(u_ode, func, t, semi::Semidiscretization)
    u = wrap_array(u_ode, semi)
    # Call `compute_coefficients` defined by the solver
    compute_coefficients!(u, func, t, mesh_equations_solver(semi)...)
end

"""
    semidiscretize(semi::Semidiscretization, tspan)

Wrap the semidiscretization `semi` as an ODE problem in the time interval `tspan`
that can be passed to `solve` from the [SciML ecosystem](https://diffeq.sciml.ai/latest/).
"""
function semidiscretize(semi::Semidiscretization, tspan)
    u0_ode = compute_coefficients(semi.initial_condition, first(tspan), semi)
    iip = true # is-inplace, i.e., we modify a vector when calling rhs!
    return ODEProblem{iip}(rhs!, u0_ode, tspan, semi)
end
