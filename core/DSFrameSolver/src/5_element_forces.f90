subroutine element_forces(dat)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: i, j, k, jj, kk
    real ::  pi, T(6,6), K_l(6,6), M_l(6,6), D_e(6), D_e_local(6), el_line(dat%acc)

    pi = 3.14159265359

    do i = 1, dat%n_elements

        call element_P2(dat, i, T, K_l, M_l)

        jj = dat%elements(i,2)
        kk = dat%elements(i,3)

        D_e = (/ dat%D_nodal(3*jj-2), dat%D_nodal(3*jj-1), dat%D_nodal(3*jj), &
                 dat%D_nodal(3*kk-2), dat%D_nodal(3*kk-1), dat%D_nodal(3*kk)  /)

        D_e_local = matmul(T,D_e)

        dat%F_all(i,:) = dat%P_ele(i,:) + matmul(K_l,D_e_local)

        call diagram_data(dat,i)
        call elastic_line(dat,i,D_e_local,el_line)

!        write (*,*) dat%el_diag

        do j = 1,dat%acc
            dat%el_diag(j,5*(i-1) + 5) = el_line(j)
        enddo

!        write (*,*) dat%el_diag

        do k = 1,dat%n_modes_to_calc
            D_e = (/ dat%eigen_v(3*jj-2,k), dat%eigen_v(3*jj-1,k), dat%eigen_v(3*jj,k), &
                     dat%eigen_v(3*kk-2,k), dat%eigen_v(3*kk-1,k), dat%eigen_v(3*kk,k)  /)

            D_e_local = matmul(T,D_e)

            call elastic_line(dat,i,D_e_local,el_line)
            do j = 1,dat%acc
                dat%eigen_v_diag((k-1)*dat%acc + j,i) = el_line(j)
            enddo
        enddo

    end do

return
end subroutine element_forces

! ==================================================================================
! ==================================================================================

subroutine find_closest(val,array,length,closest)
    implicit none

    integer :: length, i, closest
    real :: val, array(length)

    i = 1
    do while (val .GT. array(i))
        i = i + 1
    enddo

    if (abs(val - array(i)) .LE. abs(val - array(i+1))) then
        closest = i
    else
        closest = i + 1
    endif

return
end subroutine find_closest

! ==================================================================================
! ==================================================================================

subroutine diagram_data(dat,i)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: i, c, j, k, m, closest
    real ::  shear, arm, delta_x, L, val_d(3)
    real ::  xs(dat%acc), ysA(dat%acc), ysS(dat%acc), ysM(dat%acc)

    real, allocatable  :: val(:,:)
    integer, allocatable  :: idx(:)


    L = dat%El_Prop(i,8)
    delta_x = L/(dat%acc - 1)

!    write (*,*) dat%acc

    do j = 1,dat%acc
        xs(j) = (j-1)*delta_x
    enddo

    c = 0
    do k = 1, dat%n_forces_i
        if (dat%elements(i,1) .EQ. dat%forces_i(k,1)) then
            c = c + 1
        endif
    enddo

    if (allocated(idx)) deallocate (idx)
    if (allocated(val)) deallocate (val)

    allocate (idx(c))
    allocate (val(c,3))

    c = 0
    do k = 1, dat%n_forces_i
        if (dat%elements(i,1) .EQ. dat%forces_i(k,1)) then
            c = c + 1
            call find_closest(dat%forces_i(k,5)*L,xs,size(xs),closest)
            idx(c) = closest
            val(c,:) = dat%forces_i(k,2:4)
        endif
    enddo

    val_d = 0.0
    do k = 1, dat%n_forces_d
        if (dat%elements(i,1) .EQ. dat%forces_d(k,1)) then
            val_d(:) = dat%forces_d(k,2:4)
        endif
    enddo

    ysA(1) = -dat%F_all(i,1)
    ysS(1) = -dat%F_all(i,2)
    ysM(1) = -dat%F_all(i,3)
    do j = 2,dat%acc
        ysA(j) = ysA(j-1) - val_d(1)*delta_x
        ysS(j) = ysS(j-1) - val_d(2)*delta_x
        ysM(j) = ysM(j-1) - val_d(3)*delta_x
        do m = 1,size(idx)
            if (idx(m) .EQ. j) then
                ysA(j) = ysA(j) - val(m,1)
                ysS(j) = ysS(j) - val(m,2)
                ysM(j) = ysM(j) - val(m,3)
            endif
        enddo
    enddo

    shear = dat%F_all(i,2)
    arm = 0
    do j = 2,dat%acc
        arm = arm + delta_x
        ysM(j) = ysM(j) + val_d(2)*(arm**2)/2
        do m = 1,size(idx)
            if (idx(m) .LE. j) then
                ysM(j) = ysM(j) + val(m,2)*(arm-idx(m)*delta_x)
            endif
        enddo
        ysM(j) = ysM(j) + shear*arm
    enddo

    do j = 1,dat%acc
        dat%el_diag(j,5*(i-1) + 1) = xs(j)
        dat%el_diag(j,5*(i-1) + 2) = ysA(j)
        dat%el_diag(j,5*(i-1) + 3) = ysS(j)
        dat%el_diag(j,5*(i-1) + 4) = ysM(j)
    enddo

return
end subroutine diagram_data

! ==================================================================================
! ==================================================================================

subroutine elastic_line(dat,i,D_e_local,el_line)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: i, j
    real ::  delta_x, delta_y, L, D_e_local(6), A_matrix(4,4), B_vector(4), X_vector(4), ksi
    real ::  el_line(dat%acc)

    L = dat%El_Prop(i,8)
    delta_x = L/(dat%acc - 1)

    if (dat%El_Prop(i,5) .NE. 0) then
        A_matrix = reshape((/ 1. , 0. , 1.    , 0.     , &
                              0. , 1. , L     , 1.     , &
                              0. , 0. , L**2  , 2*L    , &
                              0. , 0. , L**3  , 3*L**2   /), shape(A_matrix))

        B_vector = (/ D_e_local(2), D_e_local(3), D_e_local(5), D_e_local(6)  /)

        call gauss_1(A_matrix,B_vector,X_vector,4)

        do j = 1,dat%acc
            ksi = (j-1)*delta_x
            el_line(j) = X_vector(1) + X_vector(2)*ksi + X_vector(3)*ksi**2 + X_vector(4)*ksi**3
        enddo
    else
        delta_y = (D_e_local(5) - D_e_local(2))/(dat%acc - 1)
        do j = 1,dat%acc
            el_line(j) = D_e_local(2) + (j-1)*delta_y
        enddo
    endif

return
end subroutine elastic_line

