using Trixi
using OrdinaryDiffEq

equations = ShallowWaterEquations1D(gravity_constant = 9.81, H0 = 0.8)

function initial_condition_dingemans_trixi(x, t, equations::ShallowWaterEquations1D)
    eta0 = 0.8
    A = 0.02
    # omega = 2*pi/(2.02*sqrt(2))
    k = 0.8406220896381442 # precomputed result of find_zero(k -> omega^2 - equations.gravity * k * tanh(k * eta0), 1.0) using Roots.jl
    if x[1] < -30.5 * pi / k || x[1] > -8.5 * pi / k
        h = 0.0
    else
        h = A * cos(k * x[1])
    end
    v = sqrt(equations.gravity / k * tanh(k * eta0)) * h / eta0
    if x[1] < 11.01 || x[1] >= 33.07
        b = 0.0
    elseif 11.01 <= x[1] && x[1] < 23.04
        b = 0.6 * (x[1] - 11.01) / (23.04 - 11.01)
    elseif 23.04 <= x[1] && x[1] < 27.04
        b = 0.6
    elseif 27.04 <= x[1] && x[1] < 33.07
        b = 0.6 * (33.07 - x[1]) / (33.07 - 27.04)
    else
        error("should not happen")
    end
    eta = h + eta0
    D = -b
    return Trixi.prim2cons(SVector(eta, v, b), equations)
end

initial_condition = initial_condition_dingemans_trixi

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal)
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = -138.0
coordinates_max = 46.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
tspan = (0.0, 70.0)
ode = Trixi.semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_callback = Trixi.AnalysisCallback(semi, interval = 100)

callbacks = CallbackSet(summary_callback, analysis_callback)

saveat = range(tspan..., length = 500)
sol = solve(ode, Tsit5(), reltol = 1e-7, abstol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat);
summary_callback() # print the timer summary 
