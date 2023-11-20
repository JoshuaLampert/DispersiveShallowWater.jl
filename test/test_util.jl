using Test: @test

# Use a macro to avoid world age issues when defining new initial conditions etc.
# inside an example.
"""
    @test_trixi_include(example; l2=nothing, linf=nothing, cons_error=nothing
                                change_waterheight=nothing,
                                change_velocity=nothing,
                                change_entropy=nothing,
                                change_entropy_modified=nothing,
                                lake_at_rest=nothing,
                                atol=1e-12, rtol=sqrt(eps()),
                                atol_ints=1e-11, rtol_ints=sqrt(eps()))

Test by calling `trixi_include(example; parameters...)`.
By default, only the absence of error output is checked.
If `l2`, `linf` or `cons_error` are specified, in addition the resulting L2/Linf/conservation
errors are compared approximately against these reference values, using `atol, rtol`
as absolute/relative tolerance.
If `change_waterheight`, `change_velocity`, `change_momemtum`, `change_entropy`, `change_entropy_modified`
or `lake_at_rest` are specified, in addition the resulting changes of the different errors are
compared approximately against these reference values, using `atol_ints`, `rtol_ints` as absolute/relative tolerance.
"""
macro test_trixi_include(example, args...)
    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    local cons_error = get_kwarg(args, :cons_error, nothing)
    local change_waterheight = get_kwarg(args, :change_waterheight, nothing)
    local change_velocity = get_kwarg(args, :change_velocity, nothing)
    local change_momentum = get_kwarg(args, :change_momentum, nothing)
    local change_entropy = get_kwarg(args, :change_entropy, nothing)
    local change_entropy_modified = get_kwarg(args, :change_entropy_modified, nothing)
    local lake_at_rest = get_kwarg(args, :lake_at_rest, nothing)
    local atol = get_kwarg(args, :atol, 1e-12)
    local rtol = get_kwarg(args, :rtol, sqrt(eps()))
    local atol_ints = get_kwarg(args, :atol_ints, 1e-11)
    local rtol_ints = get_kwarg(args, :rtol_ints, sqrt(eps()))

    local kwargs = Pair{Symbol, Any}[]
    for arg in args
        if (arg.head == :(=) &&
            !(arg.args[1] in (:l2, :linf, :cons_error, :change_waterheight,
                              :change_velocity, :change_momentum, :change_entropy,
                              :change_entropy_modified, :lake_at_rest,
                              :atol, :rtol, :atol_ints, :rtol_ints)))
            push!(kwargs, Pair(arg.args...))
        end
    end

    quote
        println("═"^100)
        println($example)

        # evaluate examples in the scope of the module they're called from
        @test_nowarn trixi_include(@__MODULE__, $example; $kwargs...)

        # if present, compare l2, linf and conservation errors against reference values
        if !isnothing($l2) || !isnothing($linf) || !isnothing($cons_error)
            errs = errors(analysis_callback)

            if !isnothing($l2)
                l2_measured = errs.l2_error[:, end]
                @test length($l2) == length(l2_measured)
                for (l2_expected, l2_actual) in zip($l2, l2_measured)
                    @test isapprox(l2_expected, l2_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($linf)
                linf_measured = errs.linf_error[:, end]
                @test length($linf) == length(linf_measured)
                for (linf_expected, linf_actual) in zip($linf, linf_measured)
                    @test isapprox(linf_expected, linf_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($cons_error)
                cons_error_measured = errs.conservation_error[:, end]
                @test length($cons_error) == length(cons_error_measured)
                for (conservation_error_expected, conservation_error_actual) in zip($cons_error,
                                                                                    cons_error_measured)
                    @test isapprox(conservation_error_expected, conservation_error_actual,
                                   atol = $atol, rtol = $rtol)
                end
            end
        end

        if !isnothing($change_waterheight) || !isnothing($change_velocity) ||
           !isnothing($change_momentum) ||
           !isnothing($change_entropy) || !isnothing($change_entropy_modified) ||
           !isnothing($lake_at_rest)
            ints = integrals(analysis_callback)

            if !isnothing($change_waterheight)
                waterheight_change_measured = ints.waterheight_total[end] -
                                              ints.waterheight_total[1]
                @test isapprox($change_waterheight, waterheight_change_measured,
                               atol = $atol_ints, rtol = $rtol_ints)
            end

            if !isnothing($change_velocity)
                velocity_change_measured = ints.velocity[end] - ints.velocity[1]
                @test isapprox($change_velocity, velocity_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_momentum)
                momentum_change_measured = ints.momentum[end] - ints.momentum[1]
                @test isapprox($change_momentum, momentum_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_entropy)
                entropy_change_measured = ints.entropy[end] - ints.entropy[1]
                @test isapprox($change_entropy, entropy_change_measured, atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_entropy_modified)
                entropy_modified_change_measured = ints.entropy_modified[end] -
                                                   ints.entropy_modified[1]
                @test isapprox($change_entropy_modified, entropy_modified_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($lake_at_rest)
                lake_at_rest_measured = ints.lake_at_rest_error[end]
                @test isapprox($lake_at_rest, lake_at_rest_measured, atol = $atol_ints,
                               rtol = $rtol_ints)
            end
        end
        println("═"^100)
    end
end

# Get the first value assigned to `keyword` in `args` and return `default_value`
# if there are no assignments to `keyword` in `args`.
function get_kwarg(args, keyword, default_value)
    val = default_value
    for arg in args
        if arg.head == :(=) && arg.args[1] == keyword
            val = arg.args[2]
            break
        end
    end
    return val
end

"""
    @trixi_testset "name of the testset" #= code to test #=

Similar to `@testset`, but wraps the code inside a temporary module to avoid
namespace pollution.
"""
macro trixi_testset(name, expr)
    @assert name isa String
    # TODO: `@eval` is evil
    # We would like to use
    #   mod = gensym(name)
    #   ...
    #   module $mod
    # to create new module names for every test set. However, this is not
    # compatible with the dirty hack using `@eval` to get the mapping when
    # loading structured, curvilinear meshes. Thus, we need to use a plain
    # module name here.
    quote
        local time_start = time_ns()
        @eval module TrixiTestModule
        using Test
        using DispersiveShallowWater
        include(@__FILE__)
        # We define `EXAMPLES_DIR` in (nearly) all test modules and use it to
        # get the path to the examples to be tested. However, that's not required
        # and we want to fail gracefully if it's not defined.
        try
            import ..EXAMPLES_DIR
        catch
            nothing
        end
        @testset $name $expr
        end
        local time_stop = time_ns()
        flush(stdout)
        @info("Testset "*$name*" finished in "
              *string(1.0e-9 * (time_stop - time_start))*" seconds.\n")
        nothing
    end
end
