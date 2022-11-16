mutable struct ir_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    itref_max::Cint
    acceptable_residual_relative::Float64
    acceptable_residual_absolute::Float64
    required_residual_relative::Float64
    record_residuals::Bool
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct ir_inform_type
    status::Cint
    alloc_status::Cint
    bad_alloc::NTuple{81,Cchar}
    norm_initial_residual::Float64
    norm_final_residual::Float64
end
