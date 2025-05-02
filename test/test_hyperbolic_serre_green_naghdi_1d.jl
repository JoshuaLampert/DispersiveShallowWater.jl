@testsnippet HyperbolicSerreGreenNaghdiEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d")
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007038714283663042,
                            0.006508261273448058,
                            0.0,
                            0.024517136865274798,
                            0.002141907410685252
                        ],
                        linf=[
                            0.0005088046605401519,
                            0.0036954890877776703,
                            0.0,
                            0.015022422297545818,
                            0.0013290414555349184
                        ],
                        cons_error=[
                            2.7000623958883807e-13,
                            0.00013389587974454997,
                            0.0,
                            0.005963937086921899,
                            4.502801745331908e-5
                        ],
                        change_entropy_modified=-2.3374946067633573e-7)

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl with bathymetry_mild_slope" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007038714283663042,
                            0.006508261273448058,
                            0.0,
                            0.024517136865274798,
                            0.002141907410685252
                        ],
                        linf=[
                            0.0005088046605401519,
                            0.0036954890877776703,
                            0.0,
                            0.015022422297545818,
                            0.0013290414555349184
                        ],
                        cons_error=[
                            2.7000623958883807e-13,
                            0.00013389587974454997,
                            0.0,
                            0.005963937086921899,
                            4.502801745331908e-5
                        ],
                        change_entropy_modified=-2.3374946067633573e-7)

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton_relaxation.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007041797674417557,
                            0.006510737539763134,
                            0.0,
                            0.024517447804525746,
                            0.002141928791106223
                        ],
                        linf=[
                            0.0005090662088376163,
                            0.0036987746989370907,
                            0.0,
                            0.01502004552088677,
                            0.0013289272946777064
                        ],
                        cons_error=[
                            3.126388037344441e-13,
                            0.00013409338344283483,
                            0.0,
                            0.005963706457799891,
                            4.504121848469822e-5
                        ],
                        change_entropy_modified=-5.684341886080802e-14,
                        atol=2.0e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl Jacobian" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),)

    J = @test_nowarn DispersiveShallowWater.jacobian(semi)
    @test size(J, 1) == 5 * nnodes(semi.mesh)
    @test size(J, 2) == 5 * nnodes(semi.mesh)

    # Check whether the Jacobian agrees with a finite difference approximation
    ode = semidiscretize(semi, (0.0, 0.1))
    q = ode.u0
    dq = similar(q)
    DispersiveShallowWater.rhs!(dq, q, semi, 0.0)
    h = similar(q)
    sqrt_eps = sqrt(eps())
    for i in eachindex(h)
        h[i] = sqrt_eps * (rand() - 0.5)
        q[i] += h[i]
    end
    dqh = similar(q)
    DispersiveShallowWater.rhs!(dqh, q, semi, 0.0)
    @test maximum(abs, vec(dq) + J * vec(h) - vec(dqh)) < 1.0e-11
end

@testitem "hyperbolic_serre_green_naghdi_well_balanced.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_well_balanced.jl"),
                        tspan=(0.0, 0.01),
                        l2=[
                            1.493401506906233e-14,
                            4.5636580728083475e-14,
                            4.217312212200215e-15,
                            5.399527148467818e-14,
                            4.646551952637425e-15
                        ],
                        linf=[
                            3.175237850427948e-14,
                            1.2230280990924922e-13,
                            6.661338147750939e-15,
                            1.2201879375620406e-13,
                            2.4646951146678475e-14
                        ],
                        cons_error=[
                            1.509903313490213e-14,
                            1.7256186536178874e-16,
                            2.220446049250313e-15,
                            6.788390857160712e-14,
                            2.220446049250313e-15
                        ],
                        change_entropy_modified=-3.197442310920451e-14,
                        lake_at_rest=1.833689450281284e-14)

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_manufactured.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_manufactured.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.00019628580608347483,
                            3.9074713552924354e-5,
                            0.0,
                            0.003659648461678621,
                            0.00020113904201635437
                        ],
                        linf=[
                            0.00033469307712019614,
                            8.167556996618863e-5,
                            0.0,
                            0.006470346460325516,
                            0.0003406400402869991
                        ],
                        cons_error=[0.0,
                            1.204145701494186e-5,
                            0.0,
                            1.6023987782225921,
                            4.456003127373265e-6],
                        atol=1e-7) # to make CI pass

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_dingemans.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[
                            0.22618134328164993,
                            0.7376950780256356,
                            3.389954422650516e-15,
                            0.5105085009518094,
                            0.2270170397315174
                        ],
                        linf=[
                            0.03631632774349847,
                            0.11854594481576813,
                            3.497202527569243e-15,
                            0.08108200376359903,
                            0.03645317804157178
                        ],
                        cons_error=[2.3874235921539366e-12,
                            0.0006849904339648783,
                            1.0658141036401503e-14,
                            0.0356135949814768,
                            0.00017701343821840965],
                        change_entropy=-0.0013486980525385661,
                        change_entropy_modified=-2.7995769187327824e-6)

    @test_allocations(semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_conservation.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_conservation.jl"),
                        l2=[
                            1.3601939454962908,
                            2.367812627978776,
                            3.565537895102537e-14,
                            0.8095219703329345,
                            1.3613680507825028
                        ],
                        linf=[
                            1.0010230791351136,
                            0.7870081638593119,
                            9.325873406851315e-15,
                            0.24169078965789928,
                            1.0010090485188015
                        ],
                        cons_error=[6.230038707144558e-11,
                            0.00029618251641450044,
                            4.547473508864641e-13,
                            0.00759853016721156,
                            0.0014395235644997229],
                        change_entropy=-0.20915006380346313,
                        change_entropy_modified=-0.09431070096138683,
                        atol=1e-11, # to make CI pass
                        atol_ints=4e-9) # to make CI pass

    @test_allocations(semi, sol, allocs=1_000)
end
