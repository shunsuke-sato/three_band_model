program main
  implicit none
  complex(8) :: zA(3,3),zB(3,3),lambda(3)



end program main

subroutine diag3x3(zA,zB,lambda)
  implicit none
  complex(8), intent(in) :: zA(3,3)
  complex(8), intent(out) :: zB(3,3)
  real(8), intent(out) :: lambda(3)

  real(8) :: c0,c1,c2

  c2 = -real(zA(1,1)+zA(2,2)+zA(3,3))
  c1 = real( zA(1,1)*zA(2,2) + zA(2,2)*zA(3,3) + zA(3,3)*zA(1,1) &
       -abs(zA(1,2))**2 -abs(zA(1,3))**2 -abs(zA(2,3))**2 )
  c0 = real( zA(1,1)*abs(zA(2,3))**2 + zA(2,2)*abs(zA(1,3))**2 &
       + zA(3,3)*abs(zA(1,2))**2 - zA(1,1)*zA(2,2)*zA(3,3) &
       - 2d0*conjg(zA(1,3))*zA(1,2)*zA(2,3) )

  p = c2**2 - 3d0*c1
  q = 27d0*c0/2d0 -c2**3 + 9d0*c2*c1/2d0
  t = 2d0*p**(-1.5d0)*q

end subroutine diag3x3
