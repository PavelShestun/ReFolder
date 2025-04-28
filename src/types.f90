module types
    implicit none

    type :: Atom
        character(len=4) :: name  ! CA (упрощение)
        real(8) :: x, y, z        ! Декартовы координаты
        real(8) :: mass           ! Масса атома
    end type Atom

    type :: Residue
        character(len=3) :: name  ! Тип аминокислоты (ALA, GLY, ...)
        type(Atom), allocatable :: atoms(:)  ! Массив атомов остатка
    end type Residue

    type :: ProteinStructure
        type(Residue), allocatable :: residues(:)  ! Массив остатков
        real(8), allocatable :: internal_coords(:) ! Внутренние координаты (phi, psi)
        integer :: num_residues, num_atoms
    end type ProteinStructure
end module types