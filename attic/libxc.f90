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

!> Integrate a exchange-correlation energy 
subroutine sandbox_libxc(mol,mesh)
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule
  type(tmesh), intent(in) :: mesh !< Molecular mesh

  type libxc_functional
     integer :: family ! LDA, GGA, etc.
     integer :: id     ! identifier
     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional
  type(libxc_functional) :: ifun

  integer, parameter :: funmax = 250

  character*(mline) :: line, msg, snum, word
  integer :: narg
  integer :: i, iarg
  real*8 :: e, etotal, excu, excd
  real(xc_f90_kind) :: zk
  logical :: doxc(funmax)
  integer :: idxc

  narg = command_argument_count()
  if (narg > 2) then
     doxc = .false.
     do iarg = 3, narg
        call getarg(iarg,line)
        read (line,'(I)',end=999,err=999) idxc
        doxc(idxc) = .true.
     end do
  else
     doxc = .true.
  endif

  etotal = 0d0
  do idxc = 1, funmax
     if (.not.doxc(idxc)) cycle
     ! functional 
     ifun%id = idxc
     ifun%family = xc_f90_family_from_id(ifun%id)
     if (ifun%family < 0) cycle
     if (idxc == 21) cycle ! exchange in 1D
     if (idxc == 160 .or. idxc == 182) cycle ! van Leeuwen & Baerends not implemented
     if (idxc == 207)  cycle ! BJ06 modified potential
     if (idxc == 208)  cycle ! Tran-Blaha BJ06
     if (idxc == 209)  cycle ! Rasanen, Pittali, Proetto BJ06

     call xc_f90_func_init(ifun%conf,ifun%info,ifun%id,XC_UNPOLARIZED)
     if (ifun%id == XC_LDA_C_XALPHA) then
        call xc_f90_lda_c_xalpha_set_par(ifun%conf,0.d0)
     end if
     call xc_f90_info_name(ifun%info,msg)
     line = ""
     write (snum,'(I3)') idxc
     select case (xc_f90_info_kind(ifun%info))
     case (XC_EXCHANGE)
        write (line,'("(x!",A,") ",a)') trim(adjustl(snum)), trim(adjustl(msg))
     case (XC_CORRELATION)
        write (line,'("(c!",A,") ",a)') trim(adjustl(snum)), trim(adjustl(msg))
     case (XC_EXCHANGE_CORRELATION)
        write (line,'("(xc!",A,") ",a)') trim(adjustl(snum)), trim(adjustl(msg))
     end select

     e = 0d0
     do i = 1, mesh%n
        if (mesh%rho(i,0) < 1d-30) cycle
        select case(ifun%family)
        case (XC_FAMILY_LDA)
           call xc_f90_lda_exc(ifun%conf, 1, mesh%rho(i,0), zk)
        case (XC_FAMILY_GGA)
           call xc_f90_gga_exc(ifun%conf, 1, mesh%rho(i,0), mesh%drho2(i,0), zk)
        case (XC_FAMILY_MGGA)
           call xc_f90_mgga_exc(ifun%conf, 1, mesh%rho(i,0), mesh%drho2(i,0),&
              mesh%d2rho(i,0), 0.5d0 * mesh%tau(i,0), zk)
        end select

        e = e + mesh%w(i) * mesh%rho(i,0) * zk 
     end do

     call xc_f90_func_end(ifun%conf)    
     write (iout,'(1p,E22.15,X,A)') e, trim(line)
     if (narg > 2) etotal = etotal + e
  end do
  if (narg > 2) write (iout,'(1p,E22.15,X,"Total")') etotal

  return
999 write (iout,'(/"Error in command line: postg2 file.wfx libxc [i1 i2 i3...]"/)')
  stop 1

end subroutine sandbox_libxc

