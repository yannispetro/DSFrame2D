    subroutine get_items(file_name,dat) !,nodes,elements,sections,BCs,forces)
    use datmod
    implicit none

    type(global_data)::  dat

    integer :: i
    character(len=20), intent(in) :: file_name
    character(len=6) :: val1

    val1 = 'some'

    OPEN(UNIT = 10, FILE = file_name, STATUS='UNKNOWN')
    i=0
    do while (val1 .NE. '*END')
        read (10, *) val1
        if (val1 == '*NODES') then
            dat%item(1) = i+1
        else if (val1 == '*ELEME') then
            dat%item(2) = i+1
        else if (val1 == '*MATER') then
            dat%item(3) = i+1
        else if (val1 == '*SECTI') then
            dat%item(4) = i+1
        else if (val1 == '*BC   ') then
            dat%item(5) = i+1
        else if (val1 == '*FOR_N') then
            dat%item(6) = i+1
        else if (val1 == '*FOR_I') then
            dat%item(7) = i+1
        else if (val1 == '*FOR_D') then
            dat%item(8) = i+1
        else if (val1 == '*FOR_T') then
            dat%item(9) = i+1
        else if (val1 == '*FOR_S') then
            dat%item(10) = i+1
        else if (val1 == '*FOR_F') then
            dat%item(11) = i+1
        else if (val1 == '*PARAM') then
            dat%item(12) = i+1
        end if
        i = i+1
    end do
    CLOSE(10)
    dat%item(13) = i

    return
    end subroutine get_items

    subroutine read_data(file_name,dat)
    use datmod
    implicit none

    type(global_data)::  dat

    integer :: i,j
    character(len=20), intent(in) :: file_name

    OPEN(UNIT = 10,FILE = file_name, STATUS='UNKNOWN')
    read (10, *)
    j = 1
    do i = dat%item(1)+1,dat%item(2)-1
        read (10, *) dat%nodes(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(2)+1,dat%item(3)-1
        read (10, *) dat%elements(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(3)+1,dat%item(4)-1
        read (10, *) dat%materials(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(4)+1,dat%item(5)-1
        read (10, *) dat%sections(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(5)+1,dat%item(6)-1
        read (10, *) dat%BCs(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(6)+1,dat%item(7)-1
        read (10, *) dat%forces_n(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(7)+1,dat%item(8)-1
        read (10, *) dat%forces_i(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(8)+1,dat%item(9)-1
        read (10, *) dat%forces_d(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(9)+1,dat%item(10)-1
        read (10, *) dat%forces_t(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(10)+1,dat%item(11)-1
        read (10, *) dat%forces_s(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(11)+1,dat%item(12)-1
        read (10, *) dat%forces_f(j,:)
        j = j + 1
    end do

    read (10, *)
    j = 1
    do i = dat%item(12)+1,dat%item(13)-1
        read (10, *) dat%parameters(:)
        j = j + 1
    end do
    CLOSE(10)

    return
    end subroutine read_data
