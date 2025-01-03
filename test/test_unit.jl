@testitem "Mesh1D" setup=[Setup] begin
    mesh = @test_nowarn Mesh1D(-1, 1, 10)
    @test_nowarn print(mesh)
    @test_nowarn display(mesh)
    @test ndims(mesh) == 1
    @test xmin(mesh) == -1
    @test xmax(mesh) == 1
    @test nnodes(mesh) == 10
    @test real(mesh) == Int64
end

@testitem "Solver" setup=[Setup, AdditionalImports] begin
    mesh = Mesh1D(-1.0, 1.0, 10)
    p = 3
    solver = @test_nowarn Solver(mesh, p)
    @test_nowarn print(solver)
    @test_nowarn display(solver)
    @test solver.D1 isa PeriodicDerivativeOperator
    @test solver.D2 isa PeriodicDerivativeOperator
    @test derivative_order(solver.D1) == 1
    @test derivative_order(solver.D2) == 2
    @test grid(solver) == grid(solver.D1) == grid(solver.D2)
    @test real(solver) == Float64

    D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
    uniform_mesh = UniformPeriodicMesh1D(-1.0, 1.0, 512 ÷ (p + 1))
    central = couple_discontinuously(D_legendre, uniform_mesh)
    minus = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
    plus = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
    D2 = sparse(plus) * sparse(minus)
    solver = @test_nowarn Solver(central, D2)
    @test solver.D1 isa UniformPeriodicCoupledOperator
    @test solver.D2 isa SparseMatrixCSC
    D1 = PeriodicUpwindOperators(minus, central, plus)
    solver = @test_nowarn Solver(D1, D2)
    @test solver.D1 isa PeriodicUpwindOperators
    @test solver.D2 isa SparseMatrixCSC

    D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                          accuracy_order = p, xmin = -1.0, xmax = 1.0, N = 10)
    D2 = sparse(D1.plus) * sparse(D1.minus)
    solver = @test_nowarn Solver(D1, D2)
    @test solver.D1 isa PeriodicUpwindOperators
    @test solver.D2 isa SparseMatrixCSC
    @test derivative_order(solver.D1) == 1

    p = 4 # N needs to be divisible by p
    D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
    uniform_mesh = UniformPeriodicMesh1D(-1.0, 1.0, div(12, p))
    D1 = couple_continuously(D_legendre, uniform_mesh)
    D2_legendre = legendre_second_derivative_operator(-1.0, 1.0, p + 1)
    D2 = couple_continuously(D2_legendre, uniform_mesh)
    solver = @test_nowarn Solver(D1, D2)
    @test solver.D1 isa UniformPeriodicCoupledOperator
    @test solver.D2 isa UniformPeriodicCoupledOperator
end

@testitem "Semidiscretization" setup=[Setup] begin
    equations = BBMBBMEquations1D(gravity_constant = 9.81)
    initial_condition = initial_condition_convergence_test
    boundary_conditions = boundary_condition_periodic
    mesh = Mesh1D(-1, 1, 10)
    solver = Solver(mesh, 4)
    semi = @test_nowarn Semidiscretization(mesh, equations, initial_condition, solver,
                                           boundary_conditions = boundary_conditions)
    @test_nowarn print(semi)
    @test_nowarn display(semi)
    @test ndims(semi) == ndims(mesh) == 1
    @test DispersiveShallowWater.eachnode(semi) == DispersiveShallowWater.eachnode(mesh)
    @test grid(semi) == grid(solver)
    mesh, equations, solver, cache = @test_nowarn DispersiveShallowWater.mesh_equations_solver_cache(semi)
    @test mesh == mesh
    @test equations == equations
    @test solver == solver

    equations_flat = BBMBBMEquations1D(bathymetry_type = bathymetry_flat,
                                       gravity_constant = 9.81)
    initial_condition = initial_condition_dingemans
    mesh = Mesh1D(-138, 46, 10)
    solver = Solver(mesh, 4)
    semi_flat = Semidiscretization(mesh, equations_flat, initial_condition, solver)
    @test_throws ArgumentError semidiscretize(semi_flat, (0.0, 1.0))
end

@testitem "Boundary conditions" setup=[Setup] begin
    boundary_conditions = boundary_condition_periodic
    @test_nowarn print(boundary_conditions)
    @test_nowarn display(boundary_conditions)
    boundary_conditions = boundary_condition_reflecting
    @test_nowarn print(boundary_conditions)
    @test_nowarn display(boundary_conditions)
end

@testitem "BBMEquation1D" setup=[Setup] begin
    equations = @test_nowarn @inferred BBMEquation1D(gravity_constant = 1.0)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    conversion_functions = [
        waterheight_total,
        waterheight,
        entropy,
        energy_total,
        prim2cons,
        prim2prim,
        prim2phys,
        energy_total_modified,
        entropy_modified,
        hamiltonian
    ]
    for conversion in conversion_functions
        @test DispersiveShallowWater.varnames(conversion, equations) isa Tuple
    end
    q = [42.0]
    @test @inferred(prim2prim(q, equations)) == q
    @test isapprox(@inferred(cons2prim(prim2cons(q, equations), equations)), q)
    @test @inferred(waterheight_total(q, equations)) == 42.0
    @test @inferred(waterheight(q, equations)) == 43.0
    @test @inferred(still_water_surface(q, equations)) == 0.0
    @test @inferred(prim2phys(q, equations)) == @inferred(prim2prim(q, equations))
    @testset "energy_total_modified and hamiltonian" begin
        initial_condition = initial_condition_manufactured
        boundary_conditions = boundary_condition_periodic
        mesh = @inferred Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        semi = @inferred Semidiscretization(mesh, equations, initial_condition,
                                            solver; boundary_conditions)
        q = @inferred DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                  0.0, semi)
        _, _, _, cache = @inferred DispersiveShallowWater.mesh_equations_solver_cache(semi)
        e_modified = @inferred energy_total_modified(q, equations, cache)
        e_modified_total = @inferred DispersiveShallowWater.integrate(e_modified, semi)
        @test isapprox(e_modified_total, 3.710663574870101)
        U_modified = @inferred entropy_modified(q, equations, cache)
        U_modified_total = @inferred DispersiveShallowWater.integrate(U_modified, semi)
        @test isapprox(U_modified_total, e_modified_total)
        h = @inferred hamiltonian(q, equations, cache)
        h_total = @inferred DispersiveShallowWater.integrate(h, semi)
        @test isapprox(h_total, 0.5)
    end
end

@testitem "BBMBBMEquations1D" setup=[Setup] begin
    equations = @test_nowarn @inferred BBMBBMEquations1D(gravity_constant = 9.81)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    conversion_functions = [
        waterheight_total,
        waterheight,
        velocity,
        momentum,
        discharge,
        entropy,
        energy_total,
        prim2cons,
        prim2prim,
        prim2phys,
        energy_total_modified,
        entropy_modified
    ]
    for conversion in conversion_functions
        @test DispersiveShallowWater.varnames(conversion, equations) isa Tuple
    end
    q = [42.0, 2.0, 2.0]
    @test @inferred(prim2prim(q, equations)) == q
    @test isapprox(@inferred(cons2prim(prim2cons(q, equations), equations)), q)
    @test @inferred(waterheight_total(q, equations)) == 42.0
    @test @inferred(waterheight(q, equations)) == 44.0
    @test @inferred(velocity(q, equations)) == 2.0
    @test @inferred(momentum(q, equations)) == 88.0
    @test @inferred(discharge(q, equations)) == 88.0
    @test @inferred(still_water_surface(q, equations)) == 0.0
    @test isapprox(@inferred(energy_total(q, equations)), 8740.42)
    @test @inferred(energy_total(q, equations)) == @inferred(entropy(q, equations))
    @test @inferred(prim2phys(q, equations)) == @inferred(prim2prim(q, equations))

    @testset "default implementation of energy_total_modified" begin
        initial_condition = initial_condition_convergence_test
        boundary_conditions = boundary_condition_periodic
        mesh = @inferred Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        semi = @inferred Semidiscretization(mesh, equations, initial_condition,
                                            solver; boundary_conditions)
        q = @inferred DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                  0.0, semi)
        _, _, _, cache = @inferred DispersiveShallowWater.mesh_equations_solver_cache(semi)
        e_modified = @inferred energy_total_modified(q, equations, cache)
        e_modified_total = @inferred DispersiveShallowWater.integrate(e_modified, semi)
        e_total = @inferred DispersiveShallowWater.integrate_quantity(energy_total,
                                                                      q, semi)
        @test isapprox(e_modified_total, e_total)
        U_modified = @inferred entropy_modified(q, equations, cache)
        U_modified_total = @inferred DispersiveShallowWater.integrate(U_modified, semi)
        @test isapprox(U_modified_total, e_modified_total)
    end
end

@testitem "SvaerdKalischEquations1D" setup=[Setup] begin
    equations = @test_nowarn SvaerdKalischEquations1D(gravity_constant = 9.81,
                                                      alpha = 0.0004040404040404049,
                                                      beta = 0.49292929292929294,
                                                      gamma = 0.15707070707070708)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    conversion_functions = [
        waterheight_total,
        waterheight,
        velocity,
        momentum,
        discharge,
        entropy,
        energy_total,
        prim2cons,
        prim2prim,
        prim2phys,
        energy_total_modified,
        entropy_modified
    ]
    for conversion in conversion_functions
        @test DispersiveShallowWater.varnames(conversion, equations) isa Tuple
    end
    q = [42.0, 2.0, 2.0]
    @test @inferred(prim2prim(q, equations)) == q
    @test isapprox(@inferred(cons2prim(prim2cons(q, equations), equations)), q)
    @test @inferred(waterheight_total(q, equations)) == 42.0
    @test @inferred(waterheight(q, equations)) == 44.0
    @test @inferred(velocity(q, equations)) == 2.0
    @test @inferred(momentum(q, equations)) == 88.0
    @test @inferred(discharge(q, equations)) == 88.0
    @test @inferred(still_water_surface(q, equations)) == 0.0
    @test isapprox(@inferred(energy_total(q, equations)), 8740.42)
    @test @inferred(prim2phys(q, equations)) == @inferred(prim2prim(q, equations))

    @testset "energy_total_modified" begin
        initial_condition = initial_condition_manufactured
        boundary_conditions = boundary_condition_periodic
        mesh = @inferred Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        semi = @inferred Semidiscretization(mesh, equations, initial_condition,
                                            solver; boundary_conditions)
        q = @inferred DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                  0.0, semi)
        _, _, _, cache = @inferred DispersiveShallowWater.mesh_equations_solver_cache(semi)
        e_modified = @inferred energy_total_modified(q, equations, cache)
        e_modified_total = @inferred DispersiveShallowWater.integrate(e_modified, semi)
        e_total = @inferred DispersiveShallowWater.integrate_quantity(energy_total,
                                                                      q, semi)
        @test isapprox(e_modified_total, 1450.0018635214328)
        @test isapprox(e_total, 7.405000000000001)
        U_modified = @inferred entropy_modified(q, equations, cache)
        U_modified_total = @inferred DispersiveShallowWater.integrate(U_modified, semi)
        @test isapprox(U_modified_total, e_modified_total)
    end

    @testset "reflecting boundary conditions" begin
        initial_condition = initial_condition_manufactured_reflecting
        boundary_conditions = boundary_condition_reflecting
        mesh = Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        @test_throws ArgumentError Semidiscretization(mesh, equations, initial_condition,
                                                      solver; boundary_conditions)
    end
end

@testitem "SerreGreenNaghdiEquations1D" setup=[Setup] begin
    equations = @test_nowarn @inferred SerreGreenNaghdiEquations1D(gravity_constant = 9.81)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    conversion_functions = [
        waterheight_total,
        waterheight,
        velocity,
        momentum,
        discharge,
        entropy,
        energy_total,
        prim2cons,
        prim2prim,
        prim2phys,
        energy_total_modified,
        entropy_modified
    ]
    for conversion in conversion_functions
        @test DispersiveShallowWater.varnames(conversion, equations) isa Tuple
    end
    q = [42.0, 2.0, 0.0]
    @test @inferred(prim2prim(q, equations)) == q
    @test isapprox(@inferred(cons2prim(prim2cons(q, equations), equations)), q)
    @test @inferred(waterheight_total(q, equations)) == 42.0
    @test @inferred(waterheight(q, equations)) == 42.0
    @test @inferred(velocity(q, equations)) == 2.0
    @test @inferred(momentum(q, equations)) == 84.0
    @test @inferred(discharge(q, equations)) == 84.0
    @test @inferred(still_water_surface(q, equations)) == 0.0
    @test @inferred(prim2phys(q, equations)) == @inferred(prim2prim(q, equations))

    @testset "energy_total_modified" begin
        initial_condition = initial_condition_convergence_test
        boundary_conditions = boundary_condition_periodic
        mesh = @inferred Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        semi = @inferred Semidiscretization(mesh, equations, initial_condition,
                                            solver; boundary_conditions)
        q = @inferred DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                  0.0, semi)
        _, _, _, cache = @inferred DispersiveShallowWater.mesh_equations_solver_cache(semi)
        e_modified = @inferred energy_total_modified(q, equations, cache)
        e_modified_total = @inferred DispersiveShallowWater.integrate(e_modified, semi)
        e_total = @inferred DispersiveShallowWater.integrate_quantity(energy_total,
                                                                      q, semi)
        @test isapprox(e_modified_total, 14.303587674490101)
        @test isapprox(e_total, 14.301514636021535)
        U_modified = @inferred entropy_modified(q, equations, cache)
        U_modified_total = @inferred DispersiveShallowWater.integrate(U_modified, semi)
        @test isapprox(U_modified_total, e_modified_total)
    end
end

@testitem "HyperbolicSerreGreenNaghdiEquations1D" setup=[Setup] begin
    equations = @test_nowarn @inferred HyperbolicSerreGreenNaghdiEquations1D(gravity_constant = 9.81,
                                                                             lambda = 500.0)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    conversion_functions = [
        waterheight_total,
        waterheight,
        velocity,
        momentum,
        discharge,
        entropy,
        energy_total,
        prim2cons,
        prim2prim,
        prim2phys,
        energy_total_modified,
        entropy_modified
    ]
    for conversion in conversion_functions
        @test DispersiveShallowWater.varnames(conversion, equations) isa Tuple
    end
    q = [42.0, 2.0, 0.0, -0.5, 43.0]
    @test @inferred(prim2prim(q, equations)) == q
    @test isapprox(@inferred(cons2prim(prim2cons(q, equations), equations)), q)
    @test @inferred(waterheight_total(q, equations)) == 42.0
    @test @inferred(waterheight(q, equations)) == 42.0
    @test @inferred(velocity(q, equations)) == 2.0
    @test @inferred(momentum(q, equations)) == 84.0
    @test @inferred(discharge(q, equations)) == 84.0
    @test @inferred(still_water_surface(q, equations)) == 0.0
    @test @inferred(prim2phys(q, equations)) == [42.0, 2.0, 0.0]

    @testset "energy_total_modified" begin
        initial_condition = initial_condition_soliton
        boundary_conditions = boundary_condition_periodic
        mesh = @inferred Mesh1D(-1.0, 1.0, 10)
        solver = Solver(mesh, 4)
        semi = @inferred Semidiscretization(mesh, equations, initial_condition,
                                            solver; boundary_conditions)
        q = @inferred DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                  0.0, semi)
        _, _, _, cache = @inferred DispersiveShallowWater.mesh_equations_solver_cache(semi)
        e_modified = @inferred energy_total_modified(q, equations, cache)
        e_modified_total = @inferred DispersiveShallowWater.integrate(e_modified, semi)
        e_total = @inferred DispersiveShallowWater.integrate_quantity(energy_total,
                                                                      q, semi)
        @test isapprox(e_modified_total, 14.303814990428117)
        @test isapprox(e_total, 14.301514636021535)
        U_modified = @inferred entropy_modified(q, equations, cache)
        U_modified_total = @inferred DispersiveShallowWater.integrate(U_modified, semi)
        @test isapprox(U_modified_total, e_modified_total)
    end
end

@testitem "AnalysisCallback" setup=[Setup] begin
    equations = SvaerdKalischEquations1D(gravity_constant = 9.81)
    initial_condition = initial_condition_dingemans
    boundary_conditions = boundary_condition_periodic
    mesh = Mesh1D(-1, 1, 10)
    solver = Solver(mesh, 4)
    semi = Semidiscretization(mesh, equations, initial_condition, solver,
                              boundary_conditions = boundary_conditions)
    analysis_callback = AnalysisCallback(semi; interval = 10,
                                         extra_analysis_errors = (:conservation_error,),
                                         extra_analysis_integrals = (waterheight_total,
                                                                     velocity, momentum,
                                                                     discharge, entropy,
                                                                     energy_total,
                                                                     entropy_modified,
                                                                     energy_total_modified,
                                                                     lake_at_rest_error))
    @test_nowarn print(analysis_callback)
    @test_nowarn display(analysis_callback)
end

@testitem "RelaxationCallback" setup=[Setup] begin
    relaxation_callback = RelaxationCallback(invariant = entropy)
    @test_nowarn print(relaxation_callback)
    @test_nowarn display(relaxation_callback)
end

@testitem "SummaryCallback" setup=[Setup] begin
    summary_callback = SummaryCallback()
    @test_nowarn print(summary_callback)
    @test_nowarn display(summary_callback)
end

@testitem "LinearDispersionRelation" setup=[Setup] begin
    disp_rel = LinearDispersionRelation(3.0)
    @test_nowarn print(disp_rel)
    @test_nowarn display(disp_rel)
    g = 9.81
    k = 2 * pi
    frequencies = [
        0.5660455316649682,
        0.5660455316649682,
        7.700912310929906,
        3.1189522995345467,
        3.1189522995345467
    ]
    wave_speeds = [
        0.09008894437955965,
        0.09008894437955965,
        1.2256382606017253,
        0.4963966757387569,
        0.4963966757387569
    ]

    for (i, equations) in enumerate((BBMEquation1D(gravity_constant = g),
                                     BBMBBMEquations1D(gravity_constant = g),
                                     SvärdKalischEquations1D(gravity_constant = g),
                                     SerreGreenNaghdiEquations1D(gravity_constant = g),
                                     HyperbolicSerreGreenNaghdiEquations1D(gravity_constant = g,
                                                                           lambda = 1.0)))
        @test isapprox(disp_rel(equations, k), frequencies[i])
        @test isapprox(wave_speed(disp_rel, equations, k), wave_speeds[i])
        # Add test for correct broadcasting
        @test isapprox(disp_rel.(equations, [k, k]), [frequencies[i], frequencies[i]])
        @test isapprox(wave_speed.(disp_rel, equations, [k, k]),
                       [wave_speeds[i], wave_speeds[i]])
        # For the normalized wave speed we expect c(0) = 1. Use eps() to avoid division by zero in c = omega / k
        @test isapprox(wave_speed(disp_rel, equations, eps(), normalize = true), 1.0)
    end
end

@testitem "util" setup=[Setup] begin
    @test_nowarn get_examples()

    accuracy_orders = [2, 4, 6]
    for accuracy_order in accuracy_orders
        eoc_mean_values, _ = convergence_test(default_example(), 2, N = 512,
                                              tspan = (0.0, 1.0),
                                              accuracy_order = accuracy_order)
        @test isapprox(eoc_mean_values[:l2][1], accuracy_order, atol = 0.5)
        @test isapprox(eoc_mean_values[:linf][2], accuracy_order, atol = 0.5)
        @test isapprox(eoc_mean_values[:l2][1], accuracy_order, atol = 0.5)
        @test isapprox(eoc_mean_values[:linf][2], accuracy_order, atol = 0.5)

        eoc_mean_values2, _ = convergence_test(default_example(), [512, 1024],
                                               tspan = (0.0, 1.0),
                                               accuracy_order = accuracy_order)
        for kind in (:l2, :linf), variable in (1, 2)
            eoc_mean_values[kind][variable] == eoc_mean_values2[kind][variable]
        end
    end
end
