subroutine write_output(dat)

    use datmod

    implicit none

    type(global_data)::  dat

    integer :: i

    open(unit=20, file='output.txt', ACTION="write", STATUS="replace")

    write (20, *) 'NODAL FORCES'
    write (20, *) '    Node ID   Fx               Fy               M'
    do i=1,dat%n_nodes
         write (20, *) i, dat%P_nodal(3*i-2), dat%P_nodal(3*i-1), dat%P_nodal(3*i)
    end do

    write (20, *) 'NODAL DISPLACEMENTS'
    write (20, *) '    Node ID   dx               dy               phi'
    do i=1,dat%n_nodes
         write (20, *) i, dat%D_nodal(3*i-2), dat%D_nodal(3*i-1), dat%D_nodal(3*i)
    end do

    write (20, *) 'ELEMENT FORCES'
    write (20, *) ' Element ID   Nj               Qj               Mj               Nk               Qk               Mk'
    do i=1,dat%n_elements
         write (20, *) i, dat%F_all(i,:)
    end do
    close(20)

    open(unit=30, file='diagrams.txt', ACTION="write", STATUS="replace")
    do i=1,dat%acc
         write (30, *) dat%el_diag(i,:)
    end do
    close(30)

    open(unit=100, file='eigenvectors.txt', ACTION="write", STATUS="replace")
    do i=1,dat%dmsn
         write (100, *) dat%eigen_v(i,:)
    end do
    close(100)

    open(unit=110, file='eigen_diagrams.txt', ACTION="write", STATUS="replace")
    do i=1,dat%n_modes_to_calc*dat%acc
         write (110, *) dat%eigen_v_diag(i,:)
    end do
    close(110)

end subroutine write_output
