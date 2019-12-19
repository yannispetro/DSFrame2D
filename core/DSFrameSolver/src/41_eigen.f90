subroutine eigen(dat)
    use datmod
    implicit none

    type(global_data)::  dat

    integer :: n, q, i, j, sorted
    real  :: norm, temp_l, vector(dat%dmsn)
    real, allocatable  :: K_eigen(:,:), M_eigen(:,:), e_val(:), e_vect(:,:), temp_v(:)
!
!    K_eigen = dat%K_total
!    M_eigen = dat%M_total


!    open(unit=50, file='K_M.txt', ACTION="write", STATUS="replace")
!    do i=1,n
!         write (50, *) M_eigen(i,:)
!    end do
!    close(50)

!    K_eigen = reshape((/ 2.05630374e-01,   8.89584493e-10,  -1.46171715e-06, &
!                         8.89584493e-10,   2.38374743e-02,   9.43440334e-06,  &
!                        -1.46171715e-06,   9.43440334e-06,   1.39685787e-02  /), shape(K_eigen))
!
!    M_eigen = reshape((/  0.22501692, -0.07509864, -0.05774453, &
!                         -0.07509864,  0.02569336,  0.01976284,  &
!                         -0.05774453,  0.01976284,  0.01524993  /), shape(M_eigen))

    n = dat%n_modes
!    if (dat%parameters(3) .EQ. 0) then
!        q = n
!    else
!        q = int(dat%parameters(3))
!    endif

    q = dat%n_modes_to_calc

    if (.not. allocated(K_eigen)) allocate (K_eigen(n,n))
    if (.not. allocated(M_eigen)) allocate (M_eigen(n,n))

    if (.not. allocated(e_val)) allocate (e_val(q))
    if (.not. allocated(e_vect)) allocate (e_vect(n,q))

    if (.not. allocated(temp_v)) allocate (temp_v(n))

    K_eigen = dat%eigen_K
    M_eigen = dat%eigen_M

    do i = 1,n
        do j = 1,n
            if (abs(K_eigen(i,j)) .LT. 0.01) then
                K_eigen(i,j) = 0.0
            endif
            if (abs(M_eigen(i,j)) .LT. 0.01) then
                M_eigen(i,j) = 0.0
            endif
        enddo
    enddo

    do i = 1,n-1
        do j = i+1,n
            K_eigen(j,i) = K_eigen(i,j)
            M_eigen(j,i) = M_eigen(i,j)
        enddo
    enddo

!    call positive_definite(n,K_eigen,pd)
!    write (*,*) pd

    if (dat%parameters(2) .EQ. 1) then
        call Jacobi(K_eigen, M_eigen, n, e_val, e_vect, dat%parameters)
    else
        call Subspace_Iteration(K_eigen, M_eigen, n, q, e_val, e_vect, dat%parameters)
    endif
!    call Householder_QR_InvIter(K_eigen, M_eigen, n, q, e_val, e_vect)

!    call Jacobi(K_eigen, M_eigen, n, e_val, e_vect)

!    call Subspace_Iteration(K_eigen, M_eigen, n, q, e_val, e_vect)

!
!    Normalize eigenvectors
!
    do i = 1,q
        norm = 0.0
        do j = 1,n
            if (abs(e_vect(j,i)) .GT. abs(norm)) then
                norm = abs(e_vect(j,i))
            endif
        enddo
        e_vect(:,i) = e_vect(:,i)/norm
    enddo

!
!    Sort the eigenvalues and eigenvectors from small to large eigenvalues
!
    sorted = 0
    do while (sorted .EQ. 0)
        sorted = 1
        do i = q,2,-1
            if (e_val(i) .LT. e_val(i-1)) then
                sorted = 0
                temp_l = e_val(i-1)
                e_val(i-1) = e_val(i)
                e_val(i) = temp_l
                temp_v = e_vect(:,i-1)
                e_vect(:,i-1) = e_vect(:,i)
                e_vect(:,i) = temp_v
            endif
        enddo
    enddo

!    write (*,*) e_val(1:q)
!    do i = 1,n
!        write (*,*) e_vect(i,1:q)
!    enddo

    do i = 1,q
        vector = 0.0
        vector(1:n) = e_vect(:,i)

        do j = 1,dat%dmsn
            dat%eigen_v(dat%ind_V(j),i) = vector(j)
        end do
    enddo

    deallocate (K_eigen, M_eigen, e_val, e_vect, temp_v)
!    call LDL(K_eigen,n,L,D)

    return
end subroutine eigen
