!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation
  use global_variables
  implicit none
  integer :: ikr,ikz,ik


  NKrz = NKr*(2*NKz+1)

  allocate(zCt(3,2,NKrz),eps(3,NKrz))
  allocate(kz0(-NKz:NKz),kz(-NKz:NKz),kr(NKr))
  allocate(ikr_table(NKrz),ikz_table(NKrz))
  zCt = 0d0; zCt(1,1,:) = 1d0; zCt(2,2,:) = 1d0
  eps = 0d0

! table
  ik = 0
  do ikr = 1,NKr
    do ikz = -NKz,NKz
      ik = ik + 1
      ikr_table(ik) = ikr
      ikz_table(ik) = ikz
    end do
  end do
  

  dkr = kr_max/dble(NKr)
  dkz = kz_max/dble(NKz)

  do ikz = -NKz,NKz
    kz0(ikz) = dkz*dble(ikz)
  end do
  kz = kz0

  do ikr = 1,NKr
    kr(ikr) = dkr*dble(ikr)
  end do

  do ik = 1,NKrz
    ikr = ikr_table(ik)
    ikz = ikz_table(ik)
    kz(ikz) = kz0(ikz)
    eps(1,ik) = eps_d
    eps(2,ik) = -0.5d0/mass_v*(kr(ikr)**2+kz(ikz)**2)
    eps(3,ik) = eps_g +0.5d0/mass_c*(kr(ikr)**2+kz(ikz)**2)
  end do

end subroutine preparation
