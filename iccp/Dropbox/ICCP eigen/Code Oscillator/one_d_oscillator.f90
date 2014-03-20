program one_d_oscillator

  use plot
  use plplot
  implicit none
  
  integer, parameter :: N = 400
!   integer ::
!   integer, allocatable ::
!   real(8) :: 
  real(8), allocatable :: walkers(:)
  
  allocate(walkers(N))
  
  call walk(N, walkers)
  print*, walkers
  
 contains

 subroutine walk(N, walkers)
  integer, intent(in) :: N
  integer :: i
  real(8), intent(out) :: walkers(N)
  real(8) :: r
  
  r = 0._8
  do i = 1, N
    call random_number(r)
    walkers(i) = r-0.5_8
  end do
end subroutine

end program
