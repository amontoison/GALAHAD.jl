mutable struct fit_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct fit_inform_type
    status::Cint
    alloc_status::Cint
    bad_alloc::NTuple{81,Cchar}
end
