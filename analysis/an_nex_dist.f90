!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8),parameter :: sigma = 0.1d0/27.211d0
  integer :: Nw
  real(8),allocatable :: nex_dist(:,:),nex_dist_s(:,:)
  real(8) :: Emin,Emax,dw,ww,weight
  integer :: ik,ikr,ikz,iw,iw2
  real(8) :: eps_t(3),de12

  read(*,*)Nw
  read(*,*)Emin,Emax
  dw = (Emax-Emin)/Nw
  allocate(nex_dist(0:Nw,0:3), nex_dist_s(0:Nw,0:3))

  do iw=0,Nw
    read(*,*)nex_dist(iw,:)
  end do

  nex_dist_s = 0d0

  do iw = 0,Nw
    do iw2 = 0,Nw
      weight = exp(-0.5d0*((iw-iw2)**2*(dw/sigma)**2))/sqrt(2d0*pi*sigma**2)
      nex_dist_s(iw2,:) = nex_dist_s(iw2,:) &
        + nex_dist(iw,:)*weight
    end do
  end do

  do iw = 0,Nw
    ww = Emin + iw*dw
    write(*,"(999e26.16e3)")ww,nex_dist_s(iw,:)
  end do

end program main
