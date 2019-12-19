subroutine modification_matrices(dat)

    use datmod
    implicit none

    type(global_data)::  dat

    integer :: i, j, k, cond
    real  :: theta, c_theta, s_theta, pi

    pi = 3.14159265359

    do i = 1,dat%dmsn
        dat%R_sup(i,i) = 1.0
    end do

    dat%n_uknown_P = sum(dat%BCs(:,2:4))

    k = 1
    do i = 1, dat%n_BCs

        if (dat%BCs(i,2) .EQ. 1) then
            dat%ind_V(dat%dmsn - dat%n_uknown_P + k) = 3*dat%BCs(i,1)-2
            k = k + 1
        end if
        if (dat%BCs(i,3) .EQ. 1) then
            dat%ind_V(dat%dmsn - dat%n_uknown_P + k) = 3*dat%BCs(i,1)-1
            k = k + 1
        end if
        if (dat%BCs(i,4) .EQ. 1) then
            dat%ind_V(dat%dmsn - dat%n_uknown_P + k) = 3*dat%BCs(i,1)
            k = k + 1
        end if

        if (dat%BCs(i,5) .NE. 0) then
            theta = pi*dat%BCs(i,5)/180
            c_theta = cos(theta)
            s_theta = sin(theta)

            dat%R_sup(3*dat%BCs(i,1)-2,3*dat%BCs(i,1)-2) = c_theta
            dat%R_sup(3*dat%BCs(i,1)-1,3*dat%BCs(i,1)-1) = c_theta
            dat%R_sup(3*dat%BCs(i,1)-2,3*dat%BCs(i,1)-1) = s_theta
            dat%R_sup(3*dat%BCs(i,1)-1,3*dat%BCs(i,1)-2) = -s_theta
        end if
    end do

    k = 1
    do i = 1, dat%dmsn
        cond = 0
        do j = 1, dat%n_uknown_P
            if (i .EQ. dat%ind_V(dat%dmsn-dat%n_uknown_P+j)) then
                cond = 1
            end if
        end do
        if (cond .EQ. 0) then
            dat%ind_V(k) = i
            k = k + 1
        end if
    end do

    do i = 1, dat%dmsn
        do j = 1, dat%dmsn
            if (dat%ind_V(i) .EQ. j) then
                dat%V_incident(i,j) = 1
            end if
        end do
    end do

return
end subroutine modification_matrices
