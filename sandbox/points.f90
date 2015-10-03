! -*- mode: F90 -*-
! Copyright (c) 2013-2015 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
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

!> Calculate properties at selected points
subroutine sandbox_points(mol,mesh)
  use types
  implicit none
 
  type(molecule), intent(in) :: mol !< Input molecule
  type(tmesh), intent(inout) :: mesh !< Molecular mesh

  character*(mline) :: line
  real*8 :: x(3), pro
  integer :: iprop
  type(props) :: pr
  logical :: mask(mprops)

  call getarg(3,line)
  read (line,*,end=999,err=999) iprop

  mask = .false.
  mask(iprop) = .true. 

  do while(.true.)
     read (iin,*,end=101,err=101) x
     call propts_point(mol,x,pr,mask,mesh)
     pro  = propts_select(pr,iprop)
     write (iout,'(3(F16.10,X),1p,E20.12)') x, pro
  end do

101 continue

  return
999 call error("sandbox_points","postg2 file.wfx points iprop << x0 y0 z0 ",2)

end subroutine sandbox_points
