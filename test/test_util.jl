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
    @eval module TrixiTestModule
      using Test
      using DispersiveShallowWater
      include(@__FILE__)
      @testset $name $expr
    end
    nothing
  end
end
