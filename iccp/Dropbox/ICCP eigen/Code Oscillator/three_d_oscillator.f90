program one_d_oscillator

  use plot
  use plplot
  implicit none
  
  integer, parameter :: N = 400, steps = 1000
  real(8) :: d, alpha, E
  real(8), allocatable :: walkers(3,:)
  alpha = 0.5_8
  d = 1._8
  E = 1._8
  allocate(walkers(3,N))
  
  
  call init_walkers
  call walk 
  E = sum(wavefct(walkers)**2 * localE(walkers)) / sum(wavefct(walkers)**2)
  print *, E
  
  
 contains
 
  subroutine walk
    integer :: i, j, counter
    real(8) :: kans, ran(3,N), p(N), walkers_new(N)
    
    d = 0.5_8
    ran = 0._8
    p = 0._8
    kans = 0._8
    walkers_new = 0._8
    
    do j = 1, steps
      call random_number(ran)
      ran = (ran - 0.5_8) * 2
      walkers_new = walkers + d * ran
      p = (wavefct(walkers_new) / wavefct(walkers))**2
    
    do i = 1, N
      if (p(i) .gt. 1._8) then
        walkers_new(i) = walkers(i)
        counter = counter + 1
      else
        call random_number(kans)
        if (kans .le. p(i)) then
          walkers_new(i) = walkers(i)
          counter = counter + 1
        endif
      endif
    end do
    
    if (mod (j, 100) == 0) then
      print *,counter, d
      d = d * 2 * counter/(N-1)/100
      counter = 0
    end if
    
    end do
    
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
   
 subroutine wavefct(walkers)
  real(8) :: wavefct(N), x(N)
    do i = 1, N
      r = sum(x(:,i))
      wavefct(i) = exp(-alpha*r**2)
  end function

end program
