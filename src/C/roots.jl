mutable struct roots_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    tol::Float64
    zero_coef::Float64
    zero_f::Float64
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct roots_inform_type
    status::Cint
    alloc_status::Cint
    bad_alloc::NTuple{81,Cchar}
end
