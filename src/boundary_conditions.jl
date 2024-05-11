struct BoundaryConditionPeriodic end

"""
    boundary_condition_periodic = DispersiveShallowWater.BoundaryConditionPeriodic()
A singleton struct indicating periodic boundary conditions.
"""
const boundary_condition_periodic = BoundaryConditionPeriodic()

Base.show(io::IO, ::BoundaryConditionPeriodic) = print(io, "boundary_condition_periodic")

struct BoundaryConditionReflecting end

"""
    boundary_condition_reflecting = DispersiveShallowWater.BoundaryConditionReflecting()
A singleton struct indicating reflecting boundary conditions.
"""
const boundary_condition_reflecting = BoundaryConditionReflecting()

function Base.show(io::IO, ::BoundaryConditionReflecting)
    print(io, "boundary_condition_reflecting")
end
