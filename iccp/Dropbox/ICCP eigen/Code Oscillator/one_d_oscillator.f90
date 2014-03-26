program one_d_oscillator

  implicit none
  
  integer, parameter :: N = 400, steps = 10000
  real(8) :: d, alpha, E
  real(8), allocatable :: walkers(:)
  alpha = 0.6_8
  d = 1._8
  E = 1._8
  allocate(walkers(N))
  

  call init_walkers
  call aanpassen 
  E = sum(wavefct(walkers)**2 * localE(walkers)) / sum(wavefct(walkers)**2)
  
  
 contains    
 
  subroutine aanpassen
    integer :: i, j, counter
    real(8) :: kans(N), ran(N), p(N), walkers_new(N), Etot
    
    d = 0.5_8
    ran = 0._8
    p = 0._8
    kans = 0._8
    walkers_new = 0._8
    Etot = 0._8
    
    do j = 1, steps
      call random_number(ran)
      ran = (ran - 0.5_8) * 2
      walkers_new = walkers + d * ran
      p = (wavefct(walkers_new) / wavefct(walkers))**2
    
      call random_number(kans)
      do i = 1, N
        if (p(i) .gt. kans(i)) then
          walkers(i) = walkers_new(i)
          counter = counter + 1
        end if
      end do
    
      if (mod (j, 100) == 0) then
        print *,counter, d
        d = d * 2 * counter/(N-1)/100
        counter = 0
      end if
      Etot = Etot + sum(localE(Walkers))
    end do
    print *, 'final energy', etot/steps/N
  end subroutine
  
  subroutine init_walkers
    call random_number(walkers)
    walkers = walkers - 0.5_8
  end subroutine
  
  function localE(x)
    real(8) :: localE(N), x(N)
    localE = alpha + x**2 * (0.5_8 - 2 * alpha**2)
  end function
  
  function wavefct(x)
    real(8) :: wavefct(N), x(N)
      wavefct = exp(-alpha*x**2)
  end function

end program
