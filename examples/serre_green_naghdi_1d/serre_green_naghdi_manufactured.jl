using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

###############################################################################
# Semidiscretization of the Serre-Green-Naghdi equations

# or bathymetry_variable instead of bathymetry_flat for later
equations = SerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_flat,
                                        gravity = 9.81)

initial_condition = initial_condition_manufactured
source_terms = source_terms_manufactured
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = 0.0
coordinates_max = 1.0
N = 256
mesh = Mesh1D(coordinates_min, coordinates_max, N)


# create solver with periodic upwind SBP operators of accuracy order 4
# (note that only central-type with even order are used)
D1 = upwind_operators(periodic_derivative_operator;
                      derivative_order = 1, accuracy_order = 4,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1, nothing)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions,
                          source_terms = source_terms)

###############################################################################
# Create `ODEProblem` and run the simulation

tspan = (0.0, 1.5)#(xmax(mesh) - xmin(mesh)) / sqrt(1.2 * equations.gravity)) # one period
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified))
                                                                 
callbacks = CallbackSet(analysis_callback, summary_callback)
callbacks = nothing
saveat = range(tspan..., length = 50)
sol = solve(ode, Tsit5(), abstol = 1e-10, reltol = 1e-10,
            save_everystep = false, callback = callbacks, 
           # saveat = saveat, 
            );

# plot(semi => sol)            
# sol

begin
    X = DispersiveShallowWater.grid(D1)
    i = length(sol)
    t = sol.t[i]
    h =  sol.u[i].x[1]
    u =  sol.u[i].x[2]
    p_h = plot(X, h,  title = "t = $t", xlabel = "", ylabel = "h", label="PDE solution - h")
    # plot!(p_h, X, h_a,  title = "t = $t", xlabel = "", ylabel = "h", label="Analytical solution")
    p_u = plot(X, u,  title = "", xlabel = "x", ylabel = "v", label="PDE solution - u")
    # plot!(p_u, X, u_a,  title = "", xlabel = "x", ylabel = "v", label="Analytical solution")
    plot(p_h, p_u, layout = (1, 2), size = (700, 500))
end

maximum(sol.u[end].x[1])
6.970730622952363

minimum(sol.u[end].x[2])
0.5540396435522645


anim = @animate for i in 1:1:length(sol.t)


    X = DispersiveShallowWater.grid(D1)
    t = sol.t[i]
    h =  sol.u[i].x[1]
    u =  sol.u[i].x[2]
    p_h = plot(X, h,  title = "t = $t", xlabel = "", ylabel = "h", label="PDE solution - h")
    # plot!(p_h, X, h_a,  title = "t = $t", xlabel = "", ylabel = "h", label="Analytical solution")
    p_u = plot(X, u,  title = "", xlabel = "x", ylabel = "v", label="PDE solution - u")
    # plot!(p_u, X, u_a,  title = "", xlabel = "x", ylabel = "v", label="Analytical solution")
    plot(p_h, p_u, layout = (1, 2), size = (700, 500))
    
end
gif(anim, fps = 10)
    



begin # set default plotting options for einheitliche und schöne Plots
    lw = 2; fontsize = 14; fontsize2 = 18

    default(
    grid=true, 
    box=:on,
    size = (700, 500),
    dpi = 300,
    lw = lw,
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    legendfont = font(fontsize-2,"Computer Modern"),
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    )
end