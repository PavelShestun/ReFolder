module io
    use types
    implicit none
contains
    subroutine read_pdb(filename, structure)
        character(len=*), intent(in) :: filename
        type(ProteinStructure), intent(out) :: structure
        integer :: i, io_status, num_residues
        character(len=80) :: line
        real(8) :: x, y, z

        ! Инициализация структуры
        structure%num_residues = 0
        structure%num_atoms = 0

        ! Открываем файл
        open(unit=10, file=filename, status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            print *, "Ошибка при открытии файла: ", filename
            stop
        end if

        ! Подсчитываем количество Cα-атомов
        num_residues = 0
        do
            read(10, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (line(1:4) == 'ATOM' .and. line(13:16) == ' CA ') then
                num_residues = num_residues + 1
            end if
        end do
        rewind(10)

        ! Выделяем память
        structure%num_residues = num_residues
        structure%num_atoms = num_residues  ! У нас только Cα-атомы
        allocate(structure%residues(num_residues))
        do i = 1, num_residues
            allocate(structure%residues(i)%atoms(1))  ! Только один атом (Cα) на остаток
            structure%residues(i)%atoms(1)%mass = 12.0d0  ! Масса углерода
        end do
        allocate(structure%internal_coords(2 * num_residues))

        ! Читаем координаты
        i = 0
        do
            read(10, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (line(1:4) == 'ATOM' .and. line(13:16) == ' CA ') then
                i = i + 1
                read(line(31:38), *) x
                read(line(39:46), *) y
                read(line(47:54), *) z
                structure%residues(i)%atoms(1)%x = x
                structure%residues(i)%atoms(1)%y = y
                structure%residues(i)%atoms(1)%z = z
                ! Устанавливаем начальные внутренние координаты (заглушка)
                structure%internal_coords(2*i-1) = 0.0d0  ! phi
                structure%internal_coords(2*i) = 0.0d0    ! psi
            end if
        end do

        close(10)
    end subroutine read_pdb

    subroutine write_pdb(filename, structure)
        character(len=*), intent(in) :: filename
        type(ProteinStructure), intent(in) :: structure
        integer :: i, io_status

        open(unit=11, file=filename, status='replace', action='write', iostat=io_status)
        if (io_status /= 0) then
            print *, "Ошибка при открытии файла для записи: ", filename
            stop
        end if

        do i = 1, structure%num_residues
            write(11, '(A6,I5,1X,A4,1X,A3,1X,I4,4X,3F8.3,2F6.2,10X,A1)') &
                'ATOM  ', i, 'CA  ', 'ALA', i, &
                structure%residues(i)%atoms(1)%x, &
                structure%residues(i)%atoms(1)%y, &
                structure%residues(i)%atoms(1)%z, &
                1.00, 0.00, 'C'
        end do
        write(11, '(A)') 'END'

        close(11)
    end subroutine write_pdb
end module io