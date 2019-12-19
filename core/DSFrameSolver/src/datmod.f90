module datmod
!
type :: global_data

    integer :: dmsn, item(13), n_nodes, n_elements, n_sections, n_materials, n_BCs, n_uknown_P, acc
    integer :: n_forces_n, n_forces_i, n_forces_d, n_forces_t, n_forces_s, n_forces_f
    integer :: n_modes, n_modes_to_calc

    real, allocatable    :: nodes(:,:), sections(:,:), materials(:,:)
    integer, allocatable :: elements(:,:), BCs(:,:)

    real, allocatable    :: forces_n(:,:), forces_i(:,:), forces_d(:,:)
    real, allocatable    :: forces_t(:,:), forces_s(:,:), forces_f(:,:)

    real, allocatable    :: parameters(:), El_Prop(:,:)

    real, allocatable    :: K_total(:,:), M_total(:,:)

    integer, allocatable :: V_incident(:,:), ind_V(:)
    real, allocatable    :: R_sup(:,:)

    real, allocatable    :: K_m(:,:), K_mm(:,:), M_mm(:,:)

    real, allocatable    :: D_nodal(:), P_nodal(:)

    real, allocatable    :: P_ele(:,:), F_all(:,:)

    real, allocatable    :: eigen_K(:,:), eigen_M(:,:), eigen_v(:,:)

    real, allocatable    :: el_diag(:,:), eigen_v_diag(:,:)


end type global_data

end module datmod
