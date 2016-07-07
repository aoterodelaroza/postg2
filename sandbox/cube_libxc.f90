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

!> Calculate the exchange-correlation energy density on a cube
subroutine sandbox_cubelibxc(mol,mesh)
#ifdef HAVE_LIBXC
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
#endif
  use types
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule
  type(tmesh), intent(inout) :: mesh

  character*(mline) :: line
  real*8 :: x0(3), x1(3), x(3), pro, rr, zk, sumzk, rmat(3,3)
  integer :: iprop, n(3), i, j, k, istat, narg, ix, ic, lp
  type(props) :: pr
  real*8, allocatable :: f(:,:,:)
  logical :: mask(mprops), ok
  real*8, parameter :: margin = 1d0 / .52917720859d0  ! 1 angstrom
  integer :: iff, idum

#ifdef HAVE_LIBXC
  type libxc_functional
     integer :: family ! LDA, GGA, etc.
     integer :: id     ! identifier
     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional
  type(libxc_functional) :: ifun(2)
#endif

#ifdef HAVE_LIBXC
  ! input
  call getarg(3,line)
  read (line,*,end=999,err=999) ix
  call getarg(4,line)
  read (line,*,end=999,err=999) ic
  call getarg(5,line)

  lp = 1
  ok = isreal(rr,line,lp)
  if (ok) then
     rr = rr / .52917720859
     narg = command_argument_count()
     if (narg == 11) then
        call getarg(6,line)
        read (line,*,end=999,err=999) x0(1)
        call getarg(7,line)
        read (line,*,end=999,err=999) x0(2)
        call getarg(8,line)
        read (line,*,end=999,err=999) x0(3)
        call getarg(9,line)
        read (line,*,end=999,err=999) x1(1)
        call getarg(10,line)
        read (line,*,end=999,err=999) x1(2)
        call getarg(11,line)
        read (line,*,end=999,err=999) x1(3)
        x0 = x0 / .52917720859
        x1 = x1 / .52917720859
     else
        ! limits and number of points
        x0 = mol%x(:,1)
        x1 = mol%x(:,1)
        do i = 1, mol%n
           x0 = min(x0,mol%x(:,i))
           x1 = max(x1,mol%x(:,i))
        end do
        x0 = x0 - margin
        x1 = x1 + margin
     end if
     n = int((x1-x0) / rr) + 1
     rmat = rr * eye3
  else
     open (iio,file=line,status='old')
     read (iio,*)
     read (iio,*)
     read (iio,*) idum, x0
     do i = 1, 3
        read (iio,*) n(i), rmat(:,i)
        x1 = x0 + n(i) * rmat(:,i)
     end do
     close (iio)
  endif

  allocate(f(0:n(3)-1,0:n(2)-1,0:n(1)-1),stat=istat)
  if (istat /= 0) call error('sandbox_cube','Error allocating f',2)

  ! initialize the functionals
  ifun(1)%id = ix ! exchange
  ifun(2)%id = ic ! correlation
  do iff = 1, 2
     ifun(iff)%family = xc_f90_family_from_id(ifun(iff)%id)
     call xc_f90_func_init(ifun(iff)%conf,ifun(iff)%info,ifun(iff)%id,XC_UNPOLARIZED)
  end do

  mask = .false.
  mask(1:16) = .true.

  ! calculate the points
  !$omp parallel do private(x,zk,sumzk,pr) schedule(dynamic)
  do i = 0, n(1)-1
     do j = 0, n(2)-1
        do k = 0, n(3)-1
           x = x0 + i * rmat(:,1) + j * rmat(:,2) + k * rmat(:,3)
           call propts_point(mol,x,pr,mask,mesh)
           sumzk = 0d0
           do iff = 1, 2
              if (ifun(iff)%id == 0) cycle
              select case(ifun(iff)%family)
              case (XC_FAMILY_LDA)
                 call xc_f90_lda_exc(ifun(iff)%conf, 1, pr%rho(0), zk)
                 sumzk = sumzk + zk * pr%rho(0)
              case (XC_FAMILY_GGA)
                 call xc_f90_gga_exc(ifun(iff)%conf, 1, pr%rho(0), pr%drho2(0), zk)
                 sumzk = sumzk + zk * pr%rho(0)
              case default
                 call error("sandbox_cubelibxc","only LDA or GGA",2)
              end select
           end do
           !$omp critical (writef)
           f(k,j,i) = abs(sumzk)
           !$omp end critical (writef)
        end do
     end do
  end do
  !$omp end parallel do

  do iff = 1, 2
     call xc_f90_func_end(ifun(iff)%conf)    
  end do

  ! write the cube header
  write(iout,'("title1")') 
  write(iout,'("title2")') 
  write(iout,'(I5,3(F12.6))') mol%n, x0
  do i = 1, 3
     write(iout,'(I5,3(F12.6))') n(i), rmat(:,i)
  end do
  do i = 1, mol%n
     write(iout,'(I4,F5.1,F11.6,F11.6,F11.6)') mol%z(i), 0d0, mol%x(:,i)
  end do
  do i = 0, n(1)-1
     do j = 0, n(2)-1
        write (iout,'(6(1x,e22.14))') (f(k,j,i),k=0,n(3)-1)
     enddo
  enddo

101 continue

  return

999 write (iout,'("postg2 file.wfx cube iprop step")')
  call error("sandbox_cube","wrong syntax",2)

#else
  call error("cube_libxc","postg2 not compiled with libxc",2)
#endif

end subroutine sandbox_cubelibxc
