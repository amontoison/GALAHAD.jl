mutable struct lms_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    memory_length::Cint
    method::Cint
    any_method::Bool
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct lms_time_type
    total::Float64
    setup::Float64
    form::Float64
    apply::Float64
    clock_total::Float64
    clock_setup::Float64
    clock_form::Float64
    clock_apply::Float64
end

mutable struct lms_inform_type
    status::Cint
    alloc_status::Cint
    length::Cint
    updates_skipped::Bool
    bad_alloc::NTuple{81,Cchar}
    time::lms_time_type
end
