module constants
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: size = 16               ! Lattice size (size x size)
    integer, parameter :: lsize = size - 1      !array size for lattice
    integer, parameter :: n = size*size           ! Number of spins
    real(dp), parameter :: minT = 0.5_dp          ! Minimum temperature
    real(dp), parameter :: dT = 0.1_dp            ! Temperature step
    integer, parameter :: mcs = 10000             ! Monte Carlo steps
    integer, parameter :: transient = 1000        ! Transient steps
    real(dp), parameter :: k_B = 1.0_dp           ! Boltzmann constant
end module constants

module lattice
    use constants
    implicit none
    integer :: spin(size, size)
    real(dp) :: T = 5.0_dp                        ! Current temperature
contains
    ! Initialize lattice with random spins
    subroutine initialize()
        real(dp) :: r
        integer :: i, j
        do j = 1, size
            do i = 1, size
                call random_number(r)
                spin(i,j) = merge(1, -1, r > 0.5_dp)
            end do
        end do
    end subroutine initialize

    ! Calculate energy at position (i,j)
    function energy(i, j) result(e)
        integer, intent(in) :: i, j
        integer :: left, right, up, down
        integer :: e
        
        ! Periodic boundary conditions
        left = merge(size, i-1, i == 1)
        right = merge(1, i+1, i == size)
        up = merge(1, j+1, j == size)
        down = merge(size, j-1, j == 1)
        
        e = -spin(i,j) * (spin(left,j) + spin(right,j) + spin(i,up) + spin(i,down))
    end function energy

    ! Attempt spin flip using Metropolis algorithm
    subroutine try_flip(i, j, delta_E, flipped)
        integer, intent(in) :: i, j
        real(dp), intent(out) :: delta_E
        logical, intent(out) :: flipped
        real(dp) :: r
        
        delta_E = -2.0_dp * energy(i,j)
        if (delta_E < 0.0_dp) then
            flipped = .true.
        else
            call random_number(r)
            if (r < exp(-delta_E/(k_B*T))) then
                flipped = .true.
            else
                flipped = .false.
            end if
        end if
        if (flipped) spin(i,j) = -spin(i,j)
    end subroutine try_flip

    ! Calculate total magnetization
    function magnetization() result(m)
        integer :: m
        m = sum(spin)
    end function magnetization

    ! Calculate total energy
    function total_energy() result(e)
        integer :: e
        integer :: i, j
        e = 0
        do j = 1, size
            do i = 1, size
                e = e + energy(i,j)
            end do
        end do
        e = e/2  ! Account for double counting
    end function total_energy
end module lattice

program ising_monte_carlo
    use constants
    use lattice
    implicit none
    integer :: i, j, step, mcs_count
    real(dp) :: delta_E, r
    integer :: E, M, M_abs
    logical :: flipped  ! <--- DECLARATION ADDED HERE
    real(dp) :: sum_E = 0.0_dp, sum_E2 = 0.0_dp
    real(dp) :: sum_M = 0.0_dp, sum_Mabs = 0.0_dp, sum_M2 = 0.0_dp, sum_M4 = 0.0_dp
    real(dp) :: avg_E_total, avg_E2_total, avg_M_total, avg_Mabs_total, avg_M2_total, avg_M4_total
    real(dp) :: avg_E, avg_M_val, avg_Mabs_val, avg_M2_per_spin, avg_E2_per_spin
    real(dp) :: X, X_prime, C, U_L
    integer :: seed_size, clock
    integer, allocatable :: seed(:)
    
    ! Initialize random seed
    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    call system_clock(count=clock)
    seed = clock + 37 * [(i, i=1,seed_size)]
    call random_seed(put=seed)
    
    ! Open output file
    open(unit=10, file='l16_ising_data1.dat', status='replace')
    write(10,'(a)') 'Temperature <M> <|M|> <M^2> X X'' <E> <E^2> C U_L'
    
    temperature_loop: do while (T >= minT)
        call initialize()
        
        ! Thermalization steps
        do step = 1, transient
            do mcs_count = 1, n
                call random_number(r)
                i = int(r*size) + 1
                call random_number(r)
                j = int(r*size) + 1
                call try_flip(i, j, delta_E, flipped)
            end do
        end do
        
        ! Reset sums
        sum_E = 0.0_dp; sum_E2 = 0.0_dp
        sum_M = 0.0_dp; sum_Mabs = 0.0_dp; sum_M2 = 0.0_dp; sum_M4 = 0.0_dp
        
        ! Main MC loop
        do step = 1, mcs
            do mcs_count = 1, n
                call random_number(r)
                i = int(r*size) + 1
                call random_number(r)
                j = int(r*size) + 1
                call try_flip(i, j, delta_E, flipped)
            end do
            
            ! Calculate observables
            E = total_energy()
            M = magnetization()
            M_abs = abs(M)
            
            ! Accumulate sums
            sum_E = sum_E + real(E, dp)
            sum_E2 = sum_E2 + real(E, dp)**2
            sum_M = sum_M + real(M, dp)
            sum_Mabs = sum_Mabs + real(M_abs, dp)
            sum_M2 = sum_M2 + real(M, dp)**2
            sum_M4 = sum_M4 + real(M, dp)**4
        end do
        
        ! Calculate averages of total quantities
        avg_E_total = sum_E / mcs
        avg_E2_total = sum_E2 / mcs
        avg_M_total = sum_M / mcs
        avg_Mabs_total = sum_Mabs / mcs
        avg_M2_total = sum_M2 / mcs
        avg_M4_total = sum_M4 / mcs
        
        ! Calculate per-spin quantities
        avg_E = avg_E_total / n
        avg_M_val = avg_M_total / n
        avg_Mabs_val = avg_Mabs_total / n
        avg_M2_per_spin = avg_M2_total / n**2
        avg_E2_per_spin = avg_E2_total / n**2
        
        ! Calculate derived quantities
        X = (avg_M2_total - avg_M_total**2) / (T * n)
        X_prime = (avg_M2_total - avg_Mabs_total**2) / (T * n)
        C = (avg_E2_total - avg_E_total**2) / (T**2 * n)
        U_L = 1.0_dp - avg_M4_total / (3.0_dp * avg_M2_total**2)
        
        write(10,'(10ES14.5)') T, avg_M_val, avg_Mabs_val, avg_M2_per_spin, X, X_prime, &
                               avg_E, avg_E2_per_spin, C, U_L
        
        T = T - dT
    end do temperature_loop
    
    close(10)
    print *, "Simulation complete. Data saved to ising_data.dat"
end program ising_monte_carlo
