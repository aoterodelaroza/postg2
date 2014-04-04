! -*- mode: F90 -*-

!> Write the molecular geometry to the standard output in the 
!> form of an xyz file.
subroutine sandbox_xyz(mol)
  implicit none

  type(molecule), intent(in) :: mol !< Input geometry

  integer :: i

  write (iout,*) mol%n
  write (iout,*)
  
  do i = 1, mol%n
     write (iout,'(A2,X,3(F20.10,X))') z2elem(mol%z(i)), mol%x(:,i) * 0.52917720859d0
  end do

end subroutine sandbox_xyz
