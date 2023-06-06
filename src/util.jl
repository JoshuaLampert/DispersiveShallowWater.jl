"""
    examples_dir()

Return the directory where the example files provided with DispersiveShallowWater.jl are located. If DispersiveShallowWater is
installed as a regular package (with `]add DispersiveShallowWater`), these files are read-only and should *not* be
modified. To find out which files are available, use, e.g., `readdir`.

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(DispersiveShallowWater, "examples")

# Note: We can't call the method below `DispersiveShallowWater.include` since that is created automatically
# inside `module DispersiveShallowWater` to `include` source files and evaluate them within the global scope
# of `DispersiveShallowWater`. However, users will want to evaluate in the global scope of `Main` or something
# similar to manage dependencies on their own.
"""
    trixi_include([mod::Module=Main,] example::AbstractString; kwargs...)

`include` the file `example` and evaluate its content in the global scope of module `mod`.
You can override specific assignments in `example` by supplying keyword arguments.
It's basic purpose is to make it easier to modify some parameters while running DispersiveShallowWater from the
REPL. Additionally, this is used in tests to reduce the computational burden for CI while still
providing examples with sensible default values for users.

Before replacing assignments in `example`, the keyword argument `maxiters` is inserted
into calls to `solve` and `DispersiveShallowWater.solve` with it's default value used in the SciML ecosystem
for ODEs, see https://diffeq.sciml.ai/stable/basics/common_solver_opts/#Miscellaneous.

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Examples

```jldoctest
julia> redirect_stdout(devnull) do
           trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_basic.jl"),
                         tspan=(0.0, 0.1))
           sol.t[end]
       end
0.1
```
"""
function trixi_include(mod::Module, elixir::AbstractString; kwargs...)
  Base.include(ex -> replace_assignments(insert_maxiters(ex); kwargs...), mod, elixir)
end

trixi_include(elixir::AbstractString; kwargs...) = trixi_include(Main, elixir; kwargs...)

# Helper methods used in the functions defined above, also copied from Trixi.jl

# Apply the function `f` to `expr` and all sub-expressions recursively.
walkexpr(f, expr::Expr) = f(Expr(expr.head, (walkexpr(f, arg) for arg in expr.args)...))
walkexpr(f, x) = f(x)

# Insert the keyword argument `maxiters` into calls to `solve` and `DispersiveShallowWater.solve`
# with default value `10^5` if it is not already present.
function insert_maxiters(expr)
  maxiters_default = 10^5

  expr = walkexpr(expr) do x
    if x isa Expr
      is_plain_solve = x.head === Symbol("call") && x.args[1] === Symbol("solve")
      is_trixi_solve = (x.head === Symbol("call") && x.args[1] isa Expr &&
                        x.args[1].head === Symbol(".") &&
                        x.args[1].args[1] === Symbol("DispersiveShallowWater") &&
                        x.args[1].args[2] isa QuoteNode &&
                        x.args[1].args[2].value === Symbol("solve"))

      if is_plain_solve || is_trixi_solve
        # Do nothing if `maxiters` is already set as keyword argument...
        for arg in x.args
          if arg isa Expr && arg.head === Symbol("kw") &&
             arg.args[1] === Symbol("maxiters")
            return x
          end
        end

        # ...and insert it otherwise.
        push!(x.args, Expr(Symbol("kw"), Symbol("maxiters"), maxiters_default))
      end
    end
    return x
  end

  return expr
end

# Replace assignments to `key` in `expr` by `key = val` for all `(key,val)` in `kwargs`.
function replace_assignments(expr; kwargs...)
  # replace explicit and keyword assignments
  expr = walkexpr(expr) do x
    if x isa Expr
      for (key, val) in kwargs
        if (x.head === Symbol("=") || x.head === :kw) && x.args[1] === Symbol(key)
          x.args[2] = :($val)
          # dump(x)
        end
      end
    end
    return x
  end

  return expr
end

# find a (keyword or common) assignment to `destination` in `expr`
# and return the assigned value
function find_assignment(expr, destination)
  # declare result to be able to assign to it in the closure
  local result

  # find explicit and keyword assignments
  walkexpr(expr) do x
    if x isa Expr
      if (x.head === Symbol("=") || x.head === :kw) &&
         x.args[1] === Symbol(destination)
        result = x.args[2]
        # dump(x)
      end
    end
    return x
  end

  result
end

"""
    @autoinfiltrate
    @autoinfiltrate condition::Bool

Invoke the `@infiltrate` macro of the package Infiltrator.jl to create a breakpoint for ad-hoc
interactive debugging in the REPL. If the optional argument `condition` is given, the breakpoint is
only enabled if `condition` evaluates to `true`.

As opposed to using `Infiltrator.@infiltrate` directly, this macro does not require Infiltrator.jl
to be added as a dependency to DispersiveShallowWater.jl. As a bonus, the macro will also attempt to load
the Infiltrator module if it has not yet been loaded manually.

Note: For this macro to work, the Infiltrator.jl package needs to be installed in your current Julia
environment stack.

See also: [Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl)

!!! warning "Internal use only"
    Please note that this macro is intended for internal use only. It is *not* part of the public
    API of DispersiveShallowWater.jl, and it thus can altered (or be removed) at any time without it being
    considered a breaking change.
"""
macro autoinfiltrate(condition = true)
  pkgid = Base.PkgId(Base.UUID("5903a43b-9cc3-4c30-8d17-598619ec4e9b"), "Infiltrator")
  if !haskey(Base.loaded_modules, pkgid)
    try
      Base.eval(Main, :(using Infiltrator))
    catch err
      @error "Cannot load Infiltrator.jl. Make sure it is included in your environment stack."
    end
  end
  i = get(Base.loaded_modules, pkgid, nothing)
  lnn = LineNumberNode(__source__.line, __source__.file)

  if i === nothing
    return Expr(:macrocall,
                Symbol("@warn"),
                lnn,
                "Could not load Infiltrator.")
  end

  return Expr(:macrocall,
              Expr(:., i, QuoteNode(Symbol("@infiltrate"))),
              lnn,
              esc(condition))
end