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

!> Postg2 -- a gaussian post-SCF sandbox (feb 2012).
program postg2
  use sandbox
  use meshmod
  use wfnmod
  use reader
  use types
  use param
  use tools_math
  use promolmod

  implicit none

  type(molecule) :: mol
  type(tmesh) :: mesh
  character*(mline) :: line, wfnfil
  logical :: ok
  integer :: i, idx, lp
  real*8 :: qpro, ntotal
  character*11 :: wfnext

  ! initialize
  call promol_init()
  call math_init()
  call reader_init()

  ! open the MO file
  open(unit=iimo,status='scratch',form='unformatted')
  didimo = .false.

  ! read command line
  call getarg(1,line)
  wfnfil = adjustl(trim(line))
  inquire(file=wfnfil,exist=ok)
  if (.not.ok) then
     write (iout,'("wfn file not found: ",A)') trim(wfnfil)
     stop 1
  endif

  ! read wfn and output some info
  idx = index(wfnfil,'.',.true.)
  if (wfnfil(idx+1:idx+3) == "wfx") then
     wfnext = "wfx"
     mol = readwfx(wfnfil)
  else if (wfnfil(idx+1:idx+3) == "wfn") then
     wfnext = "wfn"
     mol = readwfn(wfnfil)
  else if (wfnfil(idx+1:idx+4) == "fchk") then
     wfnext = "fck"
     mol = readfchk(wfnfil)
  else if (wfnfil(idx+1:idx+6) == "tcfchk") then
     wfnext = "tck"
     mol = readtck(wfnfil)
  else if (wfnfil(idx+1:idx+6) == "molden") then
     wfnext = "molden"
     mol = readmolden(wfnfil)
  else
     wfnext = "assumed wfn"
     mol = readwfn(wfnfil)
  endif

  ! generate the mesh, electrostatic potential
  mesh = genmesh(mol)

  ! generate the promolecular density and the in-vacuo densities
  call atomin(mol,mesh)

  ! do it, whatever it is
  call driver_sandbox(line,mol,wfnext,mesh)

  ! wrap up
  call closein()
  close(iimo)

end program postg2

