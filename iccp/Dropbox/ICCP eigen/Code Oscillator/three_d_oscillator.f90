program three_d_oscillator

  implicit none
  
  integer, parameter :: N = 400, steps = 10000
  real(8) :: d, alpha, E, test(3,5), ding(5)
  real(8), allocatable :: walkers(:,:)
  alpha = 0.6_8
  d = 1._8
  E = 1._8
  allocate(walkers(3,N))
  
  
  call init_walkers
  call metropolis 
!   E = (wavefct(walkers(:,5))**2 * localE(walkers(:,5))) / wavefct(walkers(:,5))**2)
  E = sum(wavefct(walkers)**2 * localE(walkers)) / sum(wavefct(walkers)**2)
!   print *, wavefct(walkers)**2, E
  
  
 contains
 
  subroutine metropolis
    integer :: i, j, counter, check_every
    real(8) :: kans(N), ran(3,N), p(N), walkers_new(3, N), Etot
    
    d = 0.4_8
    ran = 0._8
    p = 0._8
    kans = 0._8
    walkers_new = 0._8
    check_every = 100
    Etot = 0._8
      
    do j = 1, steps
      call random_number(ran)
      ran = (ran - 0.5_8) * 2
      walkers_new = walkers + d * ran
      p = (wavefct(walkers_new) / wavefct(walkers))**2
      
    call random_number(kans)
    do i = 1, N
      if (p(i) .gt. kans(i)) then
        walkers(:,i) = walkers_new(:,i)
        counter = counter + 1
      endif
    end do
    
    if (mod (j, check_every) == 0) then
      print *,counter/40000._8, d
      d = d * 2._8 * counter/(N-1)/check_every
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
  
!   function localE(x)
!   real(8) :: localE(1), x(3,1), rsquared(1)
!   rsquared = 0._8
!   rsquared = sum(x**2, dim=1)
!   localE = alpha + rsquared * (0.5_8 - 2._8 * alpha**2)
! end function
!    
! function wavefct(x)
!   real(8) :: x(3,1), rsquared(1), wavefct(1)
!   rsquared = 0._8
!   rsquared = sum(x**2, dim=1)
!   wavefct = exp(-alpha*rsquared)
! end function
  
  function localE(x)
    real(8) :: localE(N), x(3,N), rsquared(N)
    rsquared = 0._8
    rsquared = sum(x**2, dim=1)
    localE = alpha + rsquared * (0.5_8 - 2._8 * alpha**2)
  end function
   
  function wavefct(x)
    real(8) :: x(3,N), rsquared(N), wavefct(N)
    rsquared = 0._8
    rsquared = sum(x**2, dim=1)
    wavefct = exp(-alpha*rsquared)
  end function

end program
