struct BoundaryConditionPeriodic end

"""
    boundary_condition_periodic = DispersiveShallowWater.BoundaryConditionPeriodic()
A singleton struct indicating periodic boundary conditions.
"""
const boundary_condition_periodic = BoundaryConditionPeriodic()

Base.show(io::IO, ::BoundaryConditionPeriodic) = print(io, "boundary_condition_periodic")
