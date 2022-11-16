mutable struct sha_control_type
    f_indexing::Bool
    error::Cint
    out::Cint
    print_level::Cint
    approximation_algorithm::Cint
    dense_linear_solver::Cint
    max_sparse_degree::Cint
    space_critical::Bool
    deallocate_error_fatal::Bool
    prefix::NTuple{31,Cchar}
end

mutable struct sha_inform_type
    status::Cint
    alloc_status::Cint
    max_degree::Cint
    differences_needed::Cint
    max_reduced_degree::Cint
    bad_alloc::NTuple{81,Cchar}
end
