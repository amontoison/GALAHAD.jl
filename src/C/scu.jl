mutable struct scu_control_type
    f_indexing::Bool
end

mutable struct scu_inform_type
    alloc_status::Cint
    inertia::NTuple{3,Cint}
end
