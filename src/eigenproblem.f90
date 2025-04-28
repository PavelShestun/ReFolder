module eigenproblem
    implicit none
contains
    subroutine solve_eigenproblem_ic(hessian, kinetic, num_modes, eigenvalues, modes)
        real(8), intent(in) :: hessian(:, :), kinetic(:, :)
        integer, intent(in) :: num_modes
        real(8), allocatable, intent(out) :: eigenvalues(:), modes(:, :)
        integer :: n, info, lwork
        real(8), allocatable :: work(:), w(:)

        n = size(hessian, 1)
        allocate(eigenvalues(num_modes))
        allocate(modes(n, num_modes))
        allocate(w(n))

        ! Копируем матрицы, так как DSYGV их изменяет
        modes(:, 1:n) = hessian
        lwork = 3 * n
        allocate(work(lwork))

        ! Вызов LAPACK DSYGV
        call dsygv(1, 'V', 'U', n, modes, n, kinetic, n, w, work, lwork, info)
        if (info /= 0) then
            print *, "Ошибка в DSYGV: info =", info
            stop
        end if

        ! Выбираем первые num_modes мод
        eigenvalues = w(1:num_modes)
        modes = modes(:, 1:num_modes)
    end subroutine solve_eigenproblem_ic
end module eigenproblem