subroutine stiffness_matrix(dat)
    use datmod
    implicit none

    type(global_data)::  dat

    integer :: i, j, k
    integer :: n1, n2, sec, mat, dm
    real  :: x1, y1, x2, y2, dx, dy, L, phi, pi
    real  :: T(6,6), K_l(6,6), K_global(6,6),  M_l(6,6), M_global(6,6)

    pi = 3.14159265359
    phi = 0.0

    dm = int(dat%dmsn/dat%n_nodes)

    do i = 1,dat%n_elements

        sec = dat%elements(i,4)
        mat = int(dat%sections(sec,5))

        dat%El_Prop(i,1) = i
        dat%El_Prop(i,2) = dat%materials(mat,2) !ro !7850 # kg/m^3
        dat%El_Prop(i,3) = dat%materials(mat,3) !E  !210000000000  # N/m^2
        dat%El_Prop(i,4) = dat%sections(sec,2)  !A  !a*b # m^2
        dat%El_Prop(i,5) = dat%sections(sec,3)  !I
        dat%El_Prop(i,6) = dat%sections(sec,4)  !h
        dat%El_Prop(i,7) = dat%materials(mat,4) !a

        n1 = dat%elements(i,2)
        n2 = dat%elements(i,3)

        x1 = dat%nodes(n1,2)
        y1 = dat%nodes(n1,3)
        x2 = dat%nodes(n2,2)
        y2 = dat%nodes(n2,3)

        dx = x2 - x1
        dy = y2 - y1

        L = (dx**2 + dy**2)**0.5
        dat%El_Prop(i,8) = L

!        phi = 0.0
        if (dx .NE. 0) then
            if (dy .GE. 0 .and. dx .GT. 0) then
                phi = atan(dy/dx)
            else if (dy .GE. 0 .and. dx .LT. 0) then
                phi = pi - atan(abs(dy/dx))
            else if (dy .LE. 0 .and. dx .GT. 0) then
                phi = 2*pi - atan(abs(dy/dx))
            else if (dy .LE. 0 .and. dx .LT. 0) then
                phi = pi + atan(abs(dy/dx))
            end if
        else
            if (dy .GT. 0) then
                phi = pi/2
            else if (dy .LT. 0) then
                phi = 3*pi/2
            end if
        end if

        dat%El_Prop(i,9) = phi

        call element_P2(dat, i, T, K_l, M_l)

        K_global = matmul(matmul(transpose(T),K_l),T)
        M_global = matmul(matmul(transpose(T),M_l),T)

        do j = 1,dm
            do k = 1,dm
                dat%K_total(dm*(n1-1)+j,dm*(n1-1)+k) = dat%K_total(dm*(n1-1)+j,dm*(n1-1)+k) + K_global(j,k)
                dat%M_total(dm*(n1-1)+j,dm*(n1-1)+k) = dat%M_total(dm*(n1-1)+j,dm*(n1-1)+k) + M_global(j,k)
            end do
        end do
        do j = dm+1,2*dm
            do k = dm+1,2*dm
                dat%K_total(dm*(n2-1)+j-dm,dm*(n2-1)+k-dm) = dat%K_total(dm*(n2-1)+j-dm,dm*(n2-1)+k-dm) + K_global(j,k)
                dat%M_total(dm*(n2-1)+j-dm,dm*(n2-1)+k-dm) = dat%M_total(dm*(n2-1)+j-dm,dm*(n2-1)+k-dm) + M_global(j,k)
            end do
        end do
        do j = 1,dm
            do k = dm+1,2*dm
                dat%K_total(dm*(n1-1)+j,dm*(n2-1)+k-dm) = dat%K_total(dm*(n1-1)+j,dm*(n2-1)+k-dm) + K_global(j,k)
                dat%M_total(dm*(n1-1)+j,dm*(n2-1)+k-dm) = dat%M_total(dm*(n1-1)+j,dm*(n2-1)+k-dm) + M_global(j,k)
            end do
        end do
        do j = dm+1,2*dm
            do k = 1,dm
                dat%K_total(dm*(n2-1)+j-dm,dm*(n1-1)+k) = dat%K_total(dm*(n2-1)+j-dm,dm*(n1-1)+k) + K_global(j,k)
                dat%M_total(dm*(n2-1)+j-dm,dm*(n1-1)+k) = dat%M_total(dm*(n2-1)+j-dm,dm*(n1-1)+k) + M_global(j,k)
            end do
        end do

    end do

    return
end subroutine stiffness_matrix

subroutine element_P2(dat, el, T, K_l, M_l)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: el
    real  :: ro, E, A, Inert, phi, L, c_phi, s_phi, stA, stI, stM, T(6,6), K_l(6,6), M_l(6,6)

    ro    = dat%El_Prop(el,2)
    E     = dat%El_Prop(el,3)
    A     = dat%El_Prop(el,4)
    Inert = dat%El_Prop(el,5)
    L     = dat%El_Prop(el,8)
    phi   = dat%El_Prop(el,9)

    if (Inert .EQ. 0) then
        Inert = 0.000000000001
    endif

    c_phi = cos(phi)
    s_phi = sin(phi)
    T = reshape((/ c_phi, -s_phi, 0.0, 0.0  , 0.0   , 0.0,  &
                   s_phi, c_phi , 0.0, 0.0  , 0.0   , 0.0,  &
                   0.0  , 0.0   , 1.0, 0.0  , 0.0   , 0.0,  &
                   0.0  , 0.0   , 0.0, c_phi, -s_phi, 0.0,  &
                   0.0  , 0.0   , 0.0, s_phi, c_phi , 0.0,  &
                   0.0  , 0.0   , 0.0, 0.0  , 0.0   , 1.0  /), shape(T))

    stA = E*A/L
    stI = E*Inert/L
    K_l = reshape((/ stA , 0.0         , 0.0     , -stA, 0.0         , 0.0     ,  &
                     0.0 , 12*stI/L**2 , 6*stI/L , 0.0 , -12*stI/L**2, 6*stI/L ,  &
                     0.0 , 6*stI/L     , 4*stI   , 0.0 , -6*stI/L    , 2*stI   ,  &
                     -stA, 0.0         , 0.0     , stA , 0.0         , 0.0     ,  &
                     0.0 , -12*stI/L**2, -6*stI/L, 0.0 , 12*stI/L**2 , -6*stI/L,  &
                     0.0 , 6*stI/L     , 2*stI   , 0.0 , -6*stI/L    , 4*stI     /), shape(K_l))

    stM = ro*A*L/420.0
    M_l = reshape((/ 140.0, 0.0    , 0.0      , 70.0 , 0.0    , 0.0      ,  &
                     0.0  , 156.0  , 22.0*L   , 0.0  , 54.0   , -13.0*L  ,  &
                     0.0  , 22.0*L , 4.0*L**2 , 0.0  , 13.0*L , -3.0*L**2,  &
                     70.0 , 0.0    , 0.0      , 140.0, 0.0    , 0.0      ,  &
                     0.0  , 54.0   , 13.0*L   , 0.0  , 156.0  , -22.0*L  ,  &
                     0.0  , -13.0*L, -3.0*L**2, 0.0  , -22.0*L, 4.0*L**2   /), shape(M_l))
    M_l = M_l*stM
!    write (*,*) dat%n_elements

    return
end subroutine element_P2
