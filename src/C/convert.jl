mutable struct convert_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    transpose::Bool
    sum_duplicates::Bool
    order::Bool
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct convert_time_type
    total::Float64
    clock_total::Float64
end

mutable struct convert_inform_type
    status::Cint
    alloc_status::Cint
    duplicates::Cint
    bad_alloc::NTuple{81,Cchar}
    time::convert_time_type
end
