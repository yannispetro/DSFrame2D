subroutine solver(dat,ff)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: i, j, k, n1, n2, ss, ff, nn
    real  :: la, lb, L, phi, E, Ar, Inert, h, aT, DT, DTc
    real  :: P_ele_m(6), T(6,6), K_l(6,6), M_l(6,6), forces_s2(6), forces_s2_T(6)
    real, allocatable  :: K_ff(:,:), K_sf(:,:), K_fs(:,:), K_ss(:,:), K_ff_inv(:,:), M_ff(:,:)
    real, allocatable  :: P(:), P_m(:), P_mm(:), P_f(:), P_s(:), P_ls(:), D_s(:), D_f(:)
    real, allocatable  :: A(:,:), B(:), D_mm(:), D_m(:)


    ff = dat%dmsn - dat%n_uknown_P
    ss = dat%n_uknown_P

    if (.not. allocated(K_ff)) allocate (K_ff(ff,ff))
    if (.not. allocated(K_sf)) allocate (K_sf(ss,ff))
    if (.not. allocated(K_fs)) allocate (K_fs(ff,ss))
    if (.not. allocated(K_ss)) allocate (K_ss(ss,ss))

    if (.not. allocated(M_ff)) allocate (M_ff(ff,ff))

    if (.not. allocated(P)) allocate (P(dat%dmsn))
    if (.not. allocated(P_m)) allocate (P_m(dat%dmsn))
    if (.not. allocated(P_mm)) allocate (P_mm(dat%dmsn))

    if (.not. allocated(P_f)) allocate (P_f(ff))
    if (.not. allocated(P_ls)) allocate (P_ls(ss))

    if (.not. allocated(K_ff_inv)) allocate (K_ff_inv(ff,ff))
    if (.not. allocated(D_s)) allocate (D_s(ss))
    if (.not. allocated(D_f)) allocate (D_f(ff))
    if (.not. allocated(P_s)) allocate (P_s(ss))

    if (.not. allocated(A)) allocate (A(ff,ff))
    if (.not. allocated(B)) allocate (B(ff))

    if (.not. allocated(D_mm)) allocate (D_mm(dat%dmsn))
    if (.not. allocated(D_m)) allocate (D_m(dat%dmsn))

!    write (*,*) D_f

    do i = 1,ff
        do j = 1,ff
            K_ff(i,j) = dat%K_mm(i,j)
            M_ff(i,j) = dat%M_mm(i,j)
        end do
    end do
    do i = ff+1, dat%dmsn
        do j = 1,ff
            K_sf(i-ff,j) = dat%K_mm(i,j)
        end do
    end do
    do i = 1,ff
        do j = ff+1, dat%dmsn
            K_fs(i,j-ff) = dat%K_mm(i,j)
        end do
    end do
    do i = ff+1, dat%dmsn
        do j = ff+1, dat%dmsn
            K_ss(i-ff,j-ff) = dat%K_mm(i,j)
        end do
    end do

    P = 0.0

    do i = 1, dat%n_forces_n
        P(3*int(dat%forces_n(i,1))-2) = dat%forces_n(i,2)
        P(3*int(dat%forces_n(i,1))-1) = dat%forces_n(i,3)
        P(3*int(dat%forces_n(i,1)))   = dat%forces_n(i,4)
    end do

    dat%P_ele = 0.0

    do i = 1, dat%n_elements
        E     = dat%El_Prop(i,3)
        Ar    = dat%El_Prop(i,4)
        Inert = dat%El_Prop(i,5)
        h     = dat%El_Prop(i,6)
        aT    = dat%El_Prop(i,7)
        L     = dat%El_Prop(i,8)
        phi   = dat%El_Prop(i,9)

        do j = 1,dat%n_forces_i
            if (dat%forces_i(j,1) .EQ. i) then

                la = dat%forces_i(j,5)
                lb = 1 - dat%forces_i(j,5)

                dat%P_ele(i,1) = dat%P_ele(i,1) + (-1)*dat%forces_i(j,2)*lb
                dat%P_ele(i,2) = dat%P_ele(i,2) + (-1)*dat%forces_i(j,3)*lb**2*(3-2*lb) + 6*dat%forces_i(j,4)*la*lb/L
                dat%P_ele(i,3) = dat%P_ele(i,3) + (-1)*dat%forces_i(j,3)*lb**2*la*L + dat%forces_i(j,4)*lb*(2-3*lb)
                dat%P_ele(i,4) = dat%P_ele(i,4) + dat%forces_i(j,2)*la
                dat%P_ele(i,5) = dat%P_ele(i,5) + (-1)*dat%forces_i(j,3)*la**2*(3-2*la) + (-1)*6*dat%forces_i(j,4)*la*lb/L
                dat%P_ele(i,6) = dat%P_ele(i,6) + dat%forces_i(j,3)*la**2*lb*L + dat%forces_i(j,4)*lb*(2-3*lb)
            end if
        end do

        do j = 1,dat%n_forces_d
            if (dat%forces_d(j,1) .EQ. i) then
                dat%P_ele(i,1) = dat%P_ele(i,1) + (-1)*dat%forces_d(j,2)*L/2
                dat%P_ele(i,2) = dat%P_ele(i,2) + (-1)*dat%forces_d(j,3)*L/2
                dat%P_ele(i,3) = dat%P_ele(i,3) + (-1)*dat%forces_d(j,3)*L**2/12
                dat%P_ele(i,4) = dat%P_ele(i,4) + (-1)*dat%forces_d(j,2)*L/2
                dat%P_ele(i,5) = dat%P_ele(i,5) + (-1)*dat%forces_d(j,3)*L/2
                dat%P_ele(i,6) = dat%P_ele(i,6) + dat%forces_d(j,3)*L**2/12
            end if
        end do

        do j = 1,dat%n_forces_t
            if (dat%forces_t(j,1) .EQ. i) then
                DT = dat%forces_t(j,3) - dat%forces_t(j,2)
                DTc = (dat%forces_t(j,2) + dat%forces_t(j,3))/2 - dat%forces_t(j,4)
                dat%P_ele(i,1) = dat%P_ele(i,1) + aT*E*Ar*DTc
                dat%P_ele(i,2) = dat%P_ele(i,2) + 0.0
                dat%P_ele(i,3) = dat%P_ele(i,3) + (-1)*aT*E*Inert/h*DT
                dat%P_ele(i,4) = dat%P_ele(i,4) + (-1)*aT*E*Ar*DTc
                dat%P_ele(i,5) = dat%P_ele(i,5) + 0.0
                dat%P_ele(i,6) = dat%P_ele(i,6) + aT*E*Inert/h*DT
            end if
        end do

        call element_P2(dat, i, T, K_l, M_l)

        n1 = dat%elements(i,2)
        n2 = dat%elements(i,3)

        forces_s2 = 0.0

        do k = 1, dat%n_forces_s
            if (int(dat%forces_s(k,1)) .EQ. n1) then
                forces_s2(1) = dat%forces_s(k,2)
                forces_s2(2) = dat%forces_s(k,3)
                forces_s2(3) = dat%forces_s(k,4)
            else
                forces_s2(1) = 0.0
                forces_s2(2) = 0.0
                forces_s2(3) = 0.0
            endif

            if (int(dat%forces_s(k,1)) .EQ. n2) then
                forces_s2(4) = dat%forces_s(k,2)
                forces_s2(5) = dat%forces_s(k,3)
                forces_s2(6) = dat%forces_s(k,4)
            else
                forces_s2(4) = 0.0
                forces_s2(5) = 0.0
                forces_s2(6) = 0.0
            endif
        enddo

        forces_s2_T = matmul(T,forces_s2)
        do k = 1,size(forces_s2_T)
            if (abs(forces_s2_T(k)) .LT. 10E-9) then
                forces_s2_T(k) = 0.0
            endif
        enddo

        dat%P_ele(i,1) = dat%P_ele(i,1) + E*Ar/L*forces_s2_T(1) + (-1)*E*Ar/L*forces_s2_T(4)
        dat%P_ele(i,2) = dat%P_ele(i,2) + (-1)*12*E*Inert/L**3*(forces_s2_T(5)-forces_s2_T(2)) &
                            + 6*E*Inert/L**2*(forces_s2_T(3)+forces_s2_T(6))
        dat%P_ele(i,3) = dat%P_ele(i,3) + (-1)*6*E*Inert/L**2*(forces_s2_T(5)-forces_s2_T(2))  &
                            + 2*E*Inert/L*(2*forces_s2_T(3)+forces_s2_T(6))
        dat%P_ele(i,4) = dat%P_ele(i,4) + (-1)*E*Ar/L*forces_s2_T(1) + E*Ar/L*forces_s2_T(4)
        dat%P_ele(i,5) = dat%P_ele(i,5) + 12*E*Inert/L**3*(forces_s2_T(5)-forces_s2_T(2))      &
                            + (-1)*6*E*Inert/L**2*(forces_s2_T(3)+forces_s2_T(6))
        dat%P_ele(i,6) = dat%P_ele(i,6) + (-1)*6*E*Inert/L**2*(forces_s2_T(5)-forces_s2_T(2))  &
                            + 2*E*Inert/L*(forces_s2_T(3)+2*forces_s2_T(6))


        P_ele_m = matmul(transpose(T),dat%P_ele(i,:))

        P(3*int(dat%elements(i,2))-2) = P(3*int(dat%elements(i,2))-2) - P_ele_m(1)
        P(3*int(dat%elements(i,2))-1) = P(3*int(dat%elements(i,2))-1) - P_ele_m(2)
        P(3*int(dat%elements(i,2)))   = P(3*int(dat%elements(i,2)))   - P_ele_m(3)
        P(3*int(dat%elements(i,3))-2) = P(3*int(dat%elements(i,3))-2) - P_ele_m(4)
        P(3*int(dat%elements(i,3))-1) = P(3*int(dat%elements(i,3))-1) - P_ele_m(5)
        P(3*int(dat%elements(i,3)))   = P(3*int(dat%elements(i,3)))   - P_ele_m(6)


    end do

    if (sum(dat%BCs(:,5)) .NE. 0) then
        P_m = matmul(dat%R_sup,P)
    else
        P_m = P
    end if

    P_mm = matmul(dat%V_incident,P_m)

    P_f = P_mm(1:ff)
    P_ls = P_mm(ff+1:size(P_mm))

    !-------SOLVE---------
    D_s = 0.0

    nn = size(K_ff(1,:))

    A = K_ff
    B = P_f - matmul(K_fs,D_s)

    if (dat%parameters(1) .EQ. 1) then
        call gauss_1(A,B,D_f,nn)
!        call lu_solve(A, B)
    else
        call gauss_1(A,B,D_f,nn)
!        call lu_solve(A, B)
!        D_f = B
    endif

    P_s = matmul(K_sf,D_f) + matmul(K_ss,D_s) - P_ls
    !---------------------

    D_mm(1:ff) = D_f
    D_mm(ff+1:dat%dmsn) = D_s
    P_mm(ff+1:dat%dmsn) = P_s

    do i = 1,dat%dmsn
        D_m(dat%ind_V(i)) = D_mm(i)
        P_m(dat%ind_V(i)) = P_mm(i)
    end do

    if (sum(dat%BCs(:,5)) .NE. 0) then
        dat%D_nodal = matmul(transpose(dat%R_sup),D_m)
        dat%P_nodal = matmul(transpose(dat%R_sup),P_m)
    else
        dat%D_nodal = D_m
        dat%P_nodal = P_m
    end if

!    write (*,*) dat%P_nodal

!    write (*,*) dat%D_nodal*2.1*10**7*0.0016
!    write (*,*) P_m


    dat%eigen_K = K_ff
    dat%eigen_M = M_ff
    dat%n_modes = ff

    deallocate (K_ff, K_sf, K_fs, K_ss, K_ff_inv, M_ff)
    deallocate (P, P_m, P_mm, P_f, P_ls, P_s)
    deallocate (D_s, D_f, D_mm, D_m)
    deallocate (A, B)

return
end subroutine solver
