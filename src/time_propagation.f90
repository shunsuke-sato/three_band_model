!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it

  do it = 0,Nt
    if(myrank == 0)write(*,*)'it=',it,'/',Nt
    call dt_evolve(it)

  end do


end subroutine time_propagation
