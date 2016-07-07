! -*- mode: F90 -*-
! Copyright (c) 2013-2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
! and Axel D. Becke <axel.becke@dal.ca>
!
! postg is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
