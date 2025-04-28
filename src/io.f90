module io
    use types
    implicit none
contains
    subroutine read_pdb(filename, structure)
        character(len=*), intent(in) :: filename
        type(ProteinStructure), intent(out) :: structure
        integer :: unit, io_status, i, num_residues, num_atoms
        character(len=80) :: line
        real(8) :: x, y, z
        character(len=4) :: prev_residue_num

        ! Подсчет числа остатков и атомов
        num_residues = 0
        num_atoms = 0
        open(newunit=unit, file=filename, status='old', action='read', iostat=io_status)
        if (io_status /= 0) stop "Ошибка чтения PDB файла"
        prev_residue_num = '    '
        do
            read(unit, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (line(1:4) == 'ATOM' .and. trim(line(13:16)) == 'CA') then
                num_atoms = num_atoms + 1
                if (num_atoms == 1) then
                    num_residues = 1
                    prev_residue_num = line(23:26)
                else if (line(23:26) /= prev_residue_num) then
                    num_residues = num_residues + 1
                    prev_residue_num = line(23:26)
                end if
            end if
        end do
        close(unit)

        ! Выделение памяти
        structure%num_residues = num_residues
        structure%num_atoms = num_atoms
        allocate(structure%residues(num_residues))
        do i = 1, num_residues
            allocate(structure%residues(i)%atoms(1))  ! Только CA
            structure%residues(i)%name = 'ALA'  ! Упрощение
            structure%residues(i)%atoms(1)%name = 'CA'
            structure%residues(i)%atoms(1)%mass = 1.0d0
        end do

        ! Чтение координат
        open(newunit=unit, file=filename, status='old', action='read')
        i = 0
        do
            read(unit, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (line(1:4) == 'ATOM' .and. trim(line(13:16)) == 'CA') then
                i = i + 1
                read(line(31:38), *) x
                read(line(39:46), *) y
                read(line(47:54), *) z
                structure%residues(i)%atoms(1)%x = x
                structure%residues(i)%atoms(1)%y = y
                structure%residues(i)%atoms(1)%z = z
            end if
        end do
        close(unit)

        ! Упрощенное определение внутренних координат (заглушка)
        allocate(structure%internal_coords(2 * num_residues))
        structure%internal_coords = -60.0d0  ! Пример: phi = -60, psi = -60
    end subroutine read_pdb

    subroutine write_pdb(filename, structure)
        character(len=*), intent(in) :: filename
        type(ProteinStructure), intent(in) :: structure
        integer :: unit, i, atom_idx
        open(newunit=unit, file=filename, status='replace', action='write')
        atom_idx = 0
        do i = 1, structure%num_residues
            atom_idx = atom_idx + 1
            write(unit, '(A6,I5,1X,A4,1X,A3,1X,I4,4X,3F8.3)') &
                'ATOM  ', atom_idx, 'CA  ', 'ALA', i, &
                structure%residues(i)%atoms(1)%x, &
                structure%residues(i)%atoms(1)%y, &
                structure%residues(i)%atoms(1)%z
        end do
        write(unit, '(A)') 'END'
        close(unit)
    end subroutine write_pdb
end module io