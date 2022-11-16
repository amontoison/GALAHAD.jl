mutable struct sec_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    h_initial::Float64
    update_skip_tol::Float64
    prefix::NTuple{31,Cchar}
end

mutable struct sec_inform_type
    status::Cint
end
