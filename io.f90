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

!> Input and output routines
module io
  implicit none

  private
  public :: popwfninfo ! write info about molecule and wfn to output
  public :: read_integers ! advance lu a line and read n integers
  public :: read_reals1 ! advance lu a line and read n reals

contains

  !> Pop some information about the molecule and wavefunction to output.
  subroutine popwfninfo(mol,wfnext)
    use types
    use param
    implicit none 

    type(molecule), intent(in) :: mol !< Input molecule
    character*11, intent(in) :: wfnext !< Wavefunction file extension

    integer :: i

    write (iout,'("file          ",A)') adjustl(trim(mol%name))
    write (iout,'("filetype      ",A)') adjustl(trim(wfnext))
    write (iout,'("natoms        ",I6)') mol%n
    write (iout,'("use ECPs?     ",L5)') mol%useecp
    if (any(mol%z(:) < 0)) call error('popwfninfo','unrecognized atomic symbol',2)
    write (iout,'("# n  At           x               y               z         nr    nl")')
    do i = 1, mol%n
       write (iout,'(I4,X,A2,X,3(F15.7,X),I6,1X,I6)') &
          i, ptable(mol%z(i)), mol%x(:,i), z2nr(mol%z(i)), z2nang(mol%z(i))
    enddo
    write (iout,'("#")')
    if (mol%wfntyp == 0) then
       write (iout,'("wfn type      ","closed-shell")') 
    elseif (mol%wfntyp == 1) then
       write (iout,'("wfn type      ","open-shell")') 
    elseif (mol%wfntyp == 2) then
       write (iout,'("wfn type      ","restricted open-shell")') 
    elseif (mol%wfntyp == 3) then
       write (iout,'("wfn type      ","natural orbitals")') 
    endif
    write (iout,'("MOs           ",I6)') mol%nmo
    write (iout,'("primitives    ",I6)') mol%npri
    write (iout,'("charge        ",F6.2)') mol%charge
    write (iout,'("multiplicity  ",I6)') mol%mult

  end subroutine popwfninfo

  !> Advance LU by a single line, and read n integers from it.
  function read_integers(lu,n) result(x)
    use param
    integer, intent(in) :: lu !< LU of the file
    integer, intent(in) :: n !< Number of integers
    integer :: x(n) !< Integer list

    integer :: kk, lp, idum
    character*(mline) :: line

    kk = 0
    lp = 1
    read(lu,'(A)',end=999) line
    do while(.true.)
       if (.not.isinteger(idum,line,lp)) then
          lp = 1
          read(lu,'(A)',end=999) line
          line = adjustl(line)
          if (line(1:2) == "</") exit
       else
          kk = kk + 1
          if (kk > n) call error("read_integers","exceeded size of the array",2)
          x(kk) = idum
       endif
    enddo

    return
999 call error("read_integers","unexpected end of file",2)

  endfunction read_integers

  !> Advance LU by a single line, and read n real*8 from it.
  function read_reals1(lu,n) result(x)
    use param
    integer, intent(in) :: lu !< LU of the file
    integer, intent(in) :: n !< Number of integers
    real*8 :: x(n) !< Real*8 list

    integer :: kk, lp
    real*8 :: rdum
    character*(mline) :: line

    kk = 0
    lp = 1
    read(lu,'(A)',end=999) line
    do while(.true.)
       if (.not.isreal(rdum,line,lp)) then
          lp = 1
          read(lu,'(A)',end=999) line
          line = adjustl(line)
          if (line(1:1) == "<") exit
       else
          kk = kk + 1
          if (kk > n) call error("read_reals1","exceeded size of the array",2)
          x(kk) = rdum
       endif
    enddo

    return
999 call error("read_reals1","unexpected end of file",2)

  endfunction read_reals1

end module io
