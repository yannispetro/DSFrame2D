program main

    use datmod

    implicit none

    type(global_data)::  dat

    integer           :: ff, i
    character(len=20) :: input_file


    input_file = 'input.txt'

    call get_items(input_file,dat)

    dat%n_nodes    = dat%item(2) - dat%item(1)-1
    dat%n_elements = dat%item(3) - dat%item(2)-1
    dat%n_materials = dat%item(4) - dat%item(3)-1
    dat%n_sections = dat%item(5) - dat%item(4)-1
    dat%n_BCs      = dat%item(6) - dat%item(5)-1

    dat%n_forces_n   = dat%item(7) - dat%item(6)-1
    dat%n_forces_i   = dat%item(8) - dat%item(7)-1
    dat%n_forces_d   = dat%item(9) - dat%item(8)-1
    dat%n_forces_t   = dat%item(10) - dat%item(9)-1
    dat%n_forces_s   = dat%item(11) - dat%item(10)-1
    dat%n_forces_f   = dat%item(12) - dat%item(11)-1

    dat%dmsn = 3*dat%n_nodes

    if (.not. allocated(dat%nodes)) allocate (dat%nodes(dat%n_nodes,3))
    if (.not. allocated(dat%elements)) allocate (dat%elements(dat%n_elements,4))
    if (.not. allocated(dat%sections)) allocate (dat%sections(dat%n_sections,5))
    if (.not. allocated(dat%materials)) allocate (dat%materials(dat%n_materials,4))
    if (.not. allocated(dat%BCs)) allocate (dat%BCs(dat%n_BCs,5))

    if (dat%n_forces_n .NE. 0)then
        if (.not. allocated(dat%forces_n)) allocate (dat%forces_n(dat%n_forces_n,4))
    endif
    if (dat%n_forces_i .NE. 0)then
        if (.not. allocated(dat%forces_i)) allocate (dat%forces_i(dat%n_forces_i,5))
    endif
    if (dat%n_forces_d .NE. 0)then
        if (.not. allocated(dat%forces_d)) allocate (dat%forces_d(dat%n_forces_d,4))
    endif
    if (dat%n_forces_t .NE. 0)then
        if (.not. allocated(dat%forces_t)) allocate (dat%forces_t(dat%n_forces_t,4))
    endif
    if (dat%n_forces_s .NE. 0)then
        if (.not. allocated(dat%forces_s)) allocate (dat%forces_s(dat%n_forces_s,4))
    endif
    if (dat%n_forces_f .NE. 0)then
        if (.not. allocated(dat%forces_f)) allocate (dat%forces_f(dat%n_forces_f,4))
    endif

    if (.not. allocated(dat%parameters)) allocate (dat%parameters(12))

    if (.not. allocated(dat%El_Prop)) allocate (dat%El_Prop(dat%n_elements,9))

    call read_data(input_file, dat)

!    do i=1,dat%n_forces_n
!        write (*,*) dat%forces_n(i,:)
!    enddo

    if (.not. allocated(dat%K_total)) allocate (dat%K_total(dat%dmsn,dat%dmsn))
    if (.not. allocated(dat%M_total)) allocate (dat%M_total(dat%dmsn,dat%dmsn))
    dat%K_total = 0.0
    dat%M_total = 0.0

!    write (*,*) dat%K_total

    call stiffness_matrix(dat)

!    write (*,*) dat%K_total

!    do i = 1,dat%dmsn
!        write (*,*) dat%M_total(i,:)
!    enddo

!    write (*,*) dat%M_total

!    open(unit=500, file='K_test.txt', ACTION="write", STATUS="replace")
!    do i = 1,dat%dmsn
!        write (500,*) dat%K_total(i,:)/21/10**6/0.0016
!    enddo
!    close(500)

    if (.not. allocated(dat%R_sup)) allocate (dat%R_sup(dat%dmsn,dat%dmsn))
    if (.not. allocated(dat%ind_V)) allocate (dat%ind_V(dat%dmsn))
    if (.not. allocated(dat%V_incident)) allocate (dat%V_incident(dat%dmsn,dat%dmsn))
    dat%R_sup = 0.0
    dat%ind_V = 0
    dat%V_incident = 0

    call modification_matrices(dat)

    if (.not. allocated(dat%K_m)) allocate (dat%K_m(dat%dmsn,dat%dmsn))
    if (.not. allocated(dat%K_mm)) allocate (dat%K_mm(dat%dmsn,dat%dmsn))

    if (.not. allocated(dat%M_mm)) allocate (dat%M_mm(dat%dmsn,dat%dmsn))

    if (sum(dat%BCs(:,5)) .NE. 0) then
        dat%K_m = matmul(matmul(dat%R_sup,dat%K_total),transpose(dat%R_sup))
    else
        dat%K_m = dat%K_total
    end if

    dat%K_mm = matmul(matmul(real(dat%V_incident),dat%K_m),transpose(real(dat%V_incident)))

    dat%M_mm = matmul(matmul(real(dat%V_incident),dat%M_total),transpose(real(dat%V_incident)))

!!    do i = 1, dat%dmsn
!!        write (*,*) dat%K_mm(i,:)/21/10**6/0.0016
!!    end do
!
    if (.not. allocated(dat%D_nodal)) allocate (dat%D_nodal(dat%dmsn))
    if (.not. allocated(dat%P_nodal)) allocate (dat%P_nodal(dat%dmsn))

    if (.not. allocated(dat%P_ele)) allocate (dat%P_ele(dat%n_elements,6))

    ff = 0
    call solver(dat,ff)

    if (.not. allocated(dat%eigen_K)) allocate (dat%eigen_K(ff,ff))
    if (.not. allocated(dat%eigen_M)) allocate (dat%eigen_M(ff,ff))

    if (dat%parameters(3) .EQ. 0) then
        dat%n_modes_to_calc = dat%n_modes
    else
        dat%n_modes_to_calc = int(dat%parameters(3))
    endif

    if (.not. allocated(dat%eigen_v)) allocate (dat%eigen_v(dat%dmsn,dat%n_modes_to_calc))

    call eigen(dat)

!!    write (*,*) dat%D_nodal*2.10*10**7*0.0016
!!    write (*,*) dat%P_nodal

    if (.not. allocated(dat%F_all)) allocate (dat%F_all(dat%n_elements,6))

    dat%acc = int(dat%parameters(12) + 1)

    if (.not. allocated(dat%el_diag)) allocate (dat%el_diag(dat%acc,5*dat%n_elements))
    if (.not. allocated(dat%eigen_v_diag)) allocate (dat%eigen_v_diag(dat%acc*dat%dmsn,dat%n_elements))

    call element_forces(dat)
!!    write (*,*) dat%F_all

    call write_output(dat)

!    if (allocated(nodes, section, materials)) deallocate (nodes, section, materials)

    deallocate (dat%nodes, dat%sections, dat%materials)

    if (allocated(dat%forces_n)) deallocate (dat%forces_n)
    if (allocated(dat%forces_i)) deallocate (dat%forces_i)
    if (allocated(dat%forces_d)) deallocate (dat%forces_d)
    if (allocated(dat%forces_t)) deallocate (dat%forces_t)
    if (allocated(dat%forces_s)) deallocate (dat%forces_s)
    if (allocated(dat%forces_f)) deallocate (dat%forces_f)

    deallocate (dat%parameters)
    deallocate (dat%elements, dat%BCs, dat%V_incident, dat%ind_V)
    deallocate (dat%K_total, dat%K_m, dat%K_mm, dat%M_mm, dat%R_sup)
    deallocate (dat%D_nodal, dat%P_nodal, dat%M_total)
    deallocate (dat%El_Prop, dat%P_ele, dat%F_all, dat%el_diag)
    deallocate (dat%eigen_K, dat%eigen_M)
    deallocate (dat%eigen_v, dat%eigen_v_diag)

!    write (*,*) 'done'

    stop
end program main

