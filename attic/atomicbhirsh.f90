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

  subroutine sandbox_atomicbhirsh(mol,mesh)
    use promolmod
    use tools_math
    implicit none

    type(molecule), intent(in) :: mol
    type(tmesh), intent(in) :: mesh

    integer :: i, j
    real*8 :: ebr1, ebr2, ux, b, alf, a, x, expx
    real*8 :: etau, b1, alf1, a1, x1, r, rold, ux1
    real*8 :: dsigs, r1, r2
    integer :: ndata
    real*8, allocatable, dimension(:) :: a2, b2, c2, f2
    real*8 :: h, arho1, arho2, arho, rmid, q, r, dq, w1, w2
    real*8 :: rho, tau, dsigs, quads
    integer :: intq
    
    ! number of data points
    ndata = 200
    h = 1d0 / (ndata+1)
    rmid = 1d0

    ! spline for rho
    allocate(f2(0:ndata),a2(0:ndata),b2(0:ndata),c2(0:ndata))
    f2 = 0d0
    f2(1:ndata) = ftot(1:ndata,1)
    call spline(h,f2,a2,b2,c2,ndata,0.d0)

    r = 0d0
    rold = 0d0
    ebr1 = 0d0
    ebr2 = 0d0
    etau = 0d0
    do j = 1,2 
       do i = 1, mesh%n
          ! atom 1
          r = sqrt(dot_product(mesh%x(:,i)-mol%x(:,1),mesh%x(:,i)-mol%x(:,1)))
          q = r / (r + rmid)
          intq = int((ndata+1) * q)
          dq = q - intq * h
          arho1 = abs((f2(intq)+dq*(a2(intq)+dq*(b2(intq)+dq*c2(intq)))))/r**2
          
          ! atom 2
          r = sqrt(dot_product(mesh%x(:,i)-mol%x(:,2),mesh%x(:,i)-mol%x(:,2)))
          q = r / (r + rmid)
          intq = int((ndata+1) * q)
          dq = q - intq * h
          arho2 = abs((f2(intq)+dq*(a2(intq)+dq*(b2(intq)+dq*c2(intq)))))/r**2
          
          arho = arho1 + arho2
          w1 = arho1 / arho
          w2 = arho2 / arho
          
          ! do it
          rho = mesh%rho(i,j) * w1
          tau = mesh%tau(i,j) * w1
          dsigs = tau - 0.25d0 * pr%drho2(1) / max(pr%rho(1),1d-30)
          quads = (pr%d2rho(1) - 2d0 * dsigs) / 6d0
          call bhole(rho1,

          
          
          ! kinetic energy
          etau = etau - mesh%w(i) * 0.5d0 * mesh%tau(i,j)

          ! classic br
          b = mesh%bxdm(i,j)
          b1 = b
          alf = mesh%alf(i,j)
          alf1 = alf
          a = mesh%prefac(i,j)
          a1 = a
          x = alf * b
          x1 = x
          expx = exp(-x)
          if (b > 1d-12) then
             ux = -1d0/b * (1 - expx - 0.5d0*x*expx)
             ux1 = ux
             ebr1 = ebr1 + mesh%w(i) * 0.5d0 * mesh%rho(i,j) * ux
          endif

          ! new br 2
          r1 = sqrt(dot_product(mesh%x(:,i)-mol%x(:,1),mesh%x(:,i)-mol%x(:,1)))
          r2 = sqrt(dot_product(mesh%x(:,i)-mol%x(:,2),mesh%x(:,i)-mol%x(:,2)))
          x = 2 * r1
          expx = exp(-x)
          ux = -1d0/r1 * (1 - expx - 0.5d0*x*expx)
          x = 2 * r2
          expx = exp(-x)
          ux = ux -1d0/r2 * (1 - expx - 0.5d0*x*expx)
          ebr2 = ebr2 + mesh%w(i) * 0.5d0 * mesh%rho(i,j) * ux

          ! ! new br 1
          ! b1 = b
          ! b = mesh%bpromol(i)
          ! call bhole1(mesh%rho(i,j),b,alf,a)
          ! x = alf * b
          ! expx = exp(-x)
          ! if (b > 1d-12) then
          !    ux = -1d0/b * (1 - expx - 0.5d0*x*expx)
          !    ebr2 = ebr2 + mesh%w(i) * 0.5d0 * mesh%rho(i,j) * ux
          ! endif

          ! r1 = sqrt(dot_product(mesh%x(:,i)-mol%x(:,1),mesh%x(:,i)-mol%x(:,1)))
          ! r2 = sqrt(dot_product(mesh%x(:,i)-mol%x(:,2),mesh%x(:,i)-mol%x(:,2)))
          ! if (r1 > r2) then
          !    write (*,*) "bleh"
          !    write (*,'(1p,9(E14.6,X))') r1, r2, mesh%rho(i,j)
          !    write (*,'(1p,9(E14.6,X))') b1, b
          !    write (*,'(1p,9(E14.6,X))') a1, a
          !    write (*,'(1p,9(E14.6,X))') ux1, ux
          ! endif
       end do
    end do
    write (iout,'(1p,E22.14)') etau
    write (iout,'(1p,E22.14)') ebr1
    write (iout,'(1p,E22.14)') ebr2

  end subroutine sandbox_atomicbhirsh

