"""
    Mesh1D

Struct that holds the information for a simple homogenuous one-dimensional mesh.
"""
struct Mesh1D{RealT}
  xmin::RealT
  xmax::RealT
  N::Int

  function Mesh1D{RealT}(xmin::RealT, xmax::RealT, N::Int) where RealT
    @assert xmin < xmax
    @assert N > 0
    new(xmin, xmax, N)
  end
end

"""
    Mesh1D(xmin, xmax, N)

Create a simple homogenuous one-dimensional mesh from `xmin` to `xmax` with `N` nodes.
"""
function Mesh1D(xmin, xmax, N)
  xmin, xmax = promote(xmin, xmax)
  Mesh1D{typeof(xmin)}(xmin, xmax, N)
end

function Base.show(io::IO, mesh::Mesh1D{RealT}) where {RealT}
  print(io, "Mesh1D{", RealT, "} with length ", mesh.N)
end

function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh1D{RealT}) where {RealT}
  if get(io, :compact, false)
    show(io, mesh)
  else
    println(io, "Mesh1D{", RealT, "} ")
    println(io, "    xmin: ", mesh.xmin)
    println(io, "    xmax: ", mesh.xmax)
    print(io, "    N: ", mesh.N)
  end
end

@inline Base.ndims(mesh::Mesh1D) = 1
@inline nnodes(mesh::Mesh1D) = mesh.N
@inline Base.real(mesh::Mesh1D{RealT}) where {RealT} = RealT

xmin(mesh::Mesh1D) = mesh.xmin
xmax(mesh::Mesh1D) = mesh.xmax
