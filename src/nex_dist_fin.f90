!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine nex_dist_fin
  use global_variables
  implicit none
  integer,parameter :: Nw = 2000
  real(8) :: nex_dist(0:Nw,0:3),nex_dist_l(0:Nw,0:3)
  real(8) :: Emin,Emax,dw
  integer :: ik,ikr,ikz,iw
  real(8) :: eps_t(3),de12


  Emin = -10d0/(2d0*Ry)
  Emax =  10d0/(2d0*Ry) + eps_g
  dw = (Emax-Emin)/Nw
  nex_dist_l = 0d0
  kz = kz0

  do ik = NKrz_s,NKrz_e
    ikr = ikr_table(ik); ikz = ikz_table(ik)
    eps_t(1) = eps_d
    eps_t(2) = -0.5d0/mass_v*(kr(ikr)**2+kz(ikz)**2)
    eps_t(3) = eps_g +0.5d0/mass_c*(kr(ikr)**2+kz(ikz)**2)

! state 1
    iw = aint( (eps_t(1) - Emin)/dw )
    if(iw >= 0 .and. iw <= Nw)then
      nex_dist_l(iw,0) = nex_dist_l(iw,0) + kr(ikr)
      nex_dist_l(iw,1) = nex_dist_l(iw,1) + kr(ikr)*sum(abs(zCt(1,:,ik))**2)
    end if

! state 2
    iw = aint( (eps_t(2) - Emin)/dw )
    if(iw >= 0 .and. iw <= Nw)then
      nex_dist_l(iw,0) = nex_dist_l(iw,0) + kr(ikr)
      nex_dist_l(iw,2) = nex_dist_l(iw,2) + kr(ikr)*sum(abs(zCt(2,:,ik))**2)
    end if

! state 3
    iw = aint( (eps_t(3) - Emin)/dw )
    if(iw >= 0 .and. iw <= Nw)then
      nex_dist_l(iw,0) = nex_dist_l(iw,0) + kr(ikr)
      nex_dist_l(iw,3) = nex_dist_l(iw,3) + kr(ikr)*sum(abs(zCt(3,:,ik))**2)
    end if

  end do
  nex_dist_l=nex_dist_l*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  call MPI_ALLREDUCE(nex_dist_l,nex_dist,4*(Nw+1),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)


  if( myrank == 0)then
    open(21,file='nex_dist_fin.out')
    write(21,"(I7)")Nw
    write(21,"(2e26.16e3)")Emin,Emax
    do iw=0,Nw
      write(21,"(999e26.16e3)")nex_dist(iw,:)
    end do
    close(21)
  end if

  return
end subroutine nex_dist_fin
