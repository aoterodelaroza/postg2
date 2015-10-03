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

!> Routines that work with wavefunctions.
module wfnmod
  implicit none

  private
  public :: propts_grid ! calculate properties on the grid
  public :: propts_point ! calculate properties at a point (global and local)
  public :: mopoint ! calculate local properties at a point 
  public :: calc_exdens ! calculate the exact exchange energy density

contains

  !> Calculate properties on the grid.
  subroutine propts_grid(m,mesh,mask)
    use meshmod
    use types
    use param

    type(molecule), intent(in) :: m !< Molecule
    type(tmesh), intent(inout) :: mesh !< Mesh, with all the properties set.
    logical, intent(in) :: mask(mprops) !< Properties mask: only calculate true.

    integer :: istat, i, j
    real*8 :: prho(2), drho2(2), d2rho(2), taup(2), pb(2), phi(m%nmo)
    type(props) :: pr
    logical :: lmask(mprops)
    real*8 :: wphi(mesh%n,m%nmo), quads, dsigs, uxps

    ! electrostatic potential -> needs the total density
    lmask = mask
    if (mask(20)) lmask(1) = .true.

    ! exchange energy density -> needs the total density
    lmask = mask
    if (mask(21)) lmask(1) = .true.

    ! inverse BR hole normalization: needs the exdens and all the stuff
    lmask = mask
    if (mask(22)) then
       lmask(1:16) = .true.
       lmask(21) = .true.
    end if

    ! only calculate unknown properties
    lmask = lmask .and..not.mesh%isthere

    ! allocate
    call mesh_allocate(mesh,lmask)

    ! calculate the properties from the orbitals
    if (any(lmask(1:19))) then
       !$omp parallel do private(pr,phi) schedule(dynamic)
       do i = 1, mesh%n
          call mopoint(m,mesh%x(:,i),pr,lmask,phi)
          !$omp critical (write)
          if (any(lmask(1:4))) mesh%rho(i,:) = pr%rho
          if (any(lmask(5:8))) mesh%drho2(i,:) = pr%drho2
          if (any(lmask(9:12))) mesh%d2rho(i,:) = pr%d2rho
          if (any(lmask(13:16))) mesh%tau(i,:) = pr%taup
          if (lmask(17)) mesh%bxdm(i,:) = pr%bxdm
          if (lmask(18)) mesh%alf(i,:) = pr%alf
          if (lmask(19)) mesh%prefac(i,:) = pr%prefac
          wphi(i,:) = phi
          !$omp end critical (write)
       enddo ! i = 1, mesh%n
       !$omp end parallel do

       ! write the imo file
       rewind(iimo)
       do i = 1, m%nmo
          write(iimo) wphi(:,i)
       end do
       didimo = .true.
    endif

    ! calculate the electrostatic potential on the grid
    if (lmask(20)) then
       call poiss_grid(m,mesh,0)
    end if

    ! calculate the exchange energy density
    if (lmask(21)) then
       call calc_exdens(m,mesh)
    endif

    if (lmask(22)) then
       call calc_xlns(mesh)
    endif

    ! update what has been calculated
    mesh%isthere = mesh%isthere .or. lmask

  end subroutine propts_grid

  !> Calculate properties at a point, local and global properties. 
  subroutine propts_point(m,x,pr,mask,mesh)
    use tools, only: bhole
    use types
    use param
    implicit none

    type(molecule), intent(in) :: m !< Molecule
    real*8, intent(in) :: x(3) !< Point
    type(props), intent(out) :: pr !< Output properties
    logical, intent(in) :: mask(mprops)
    type(tmesh), intent(inout), optional :: mesh

    logical :: lmask(mprops)

    if (any(mask(1:19))) then
       ! basic properties
       call mopoint(m,x,pr,mask)

    else if (mask(20)) then
       ! electrostatic potential
       if (.not.present(mesh)) call error("propts_point","vel requires mesh",2)

       ! activate electrostatic potential and electron density
       lmask = mask_none()
       lmask(1) = .true.
       lmask(20) = .true.

       ! make sure all we need is allocated
       call mesh_allocate(mesh,lmask)

       ! calculate density and electrostatic potential on the grid
       call propts_grid(m,mesh,lmask)

       call error('propts_point','vel at a point not implemented',2)
    end if

  end subroutine propts_point

  !> Calculate basic properties on a point. Only local properties.
  subroutine mopoint(m,x,pr,mask,ophi)
    use tools, only: bhole
    use types
    use param
    implicit none
    
    type(molecule), intent(in) :: m !< Molecule
    real*8, intent(in) :: x(3) !< Point
    type(props), intent(out) :: pr !< Output properties
    logical, intent(in) :: mask(mprops)
    real*8, intent(out), optional :: ophi(m%nmo)

    integer :: iat
    integer :: i, j, nn, ityp, ipri, ipria, ix, l(3)
    integer :: imo, nspin, n0(2), n1(2), nmo1, istat
    real*8 :: al, x0(3), ex, xl(3,0:2), xl2
    real*8 :: phi(m%nmo,10), gg(3), hh(3), quads
    real*8 :: dsigs, aocc
    logical :: ldopri(m%npri,10)
    real*8 :: chi(m%npri,10), maxc(m%npri), dd(3,m%n), d2(m%n)

    real*8, parameter :: cutoff_pri = 1d-15
    integer, parameter :: li(3,56) = reshape((/&
      0,0,0, & ! s
      1,0,0, 0,1,0, 0,0,1, & ! p
      2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1, & !d
      3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, &
      1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f
      4,0,0, 0,4,0, 0,0,4, 3,1,0, 3,0,1, 1,3,0, 0,3,1, 1,0,3,&
      0,1,3, 2,2,0, 2,0,2, 0,2,2, 2,1,1, 1,2,1, 1,1,2,& ! g
      0,0,5, 0,1,4, 0,2,3, 0,3,2, 0,4,1, 0,5,0, 1,0,4, 1,1,3,&
      1,2,2, 1,3,1, 1,4,0, 2,0,3, 2,1,2, 2,2,1, 2,3,0, 3,0,2,&
      3,1,1, 3,2,0, 4,0,1, 4,1,0, 5,0,0/),shape(li)) ! h

    ! identify the max coefficients
    maxc = 0d0
    do imo = 1, m%nmo
       do ipri = 1, m%npri
          maxc(ipri) = max(maxc(ipri),abs(m%c(imo,ipri)))
       enddo
    enddo

    ! calculate distances
    do iat = 1, m%n
       dd(:,iat) = x - m%x(:,iat)
       d2(iat) = dd(1,iat)*dd(1,iat)+dd(2,iat)*dd(2,iat)+dd(3,iat)*dd(3,iat)
    enddo

    do ipri = 1, m%npri
       ityp = m%itype(ipri)
       iat = m%icenter(ipri)
       al = m%e(ipri)
       ex = exp(-al * d2(iat))

       l = li(1:3,ityp)
       do ix = 1, 3
          if (l(ix) == 0) then
             xl(ix,0) = 1d0
             xl(ix,1) = 0d0
             xl(ix,2) = 0d0
          else if (l(ix) == 1) then
             xl(ix,0) = dd(ix,iat)
             xl(ix,1) = 1d0
             xl(ix,2) = 0d0
          else if (l(ix) == 2) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 2d0 * dd(ix,iat)
             xl(ix,2) = 2d0
          else if (l(ix) == 3) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
             xl(ix,2) = 6d0 * dd(ix,iat)
          else if (l(ix) == 4) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2
             xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
             xl(ix,2) = 12d0 * xl2
          else if (l(ix) == 5) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2 * dd(ix,iat)
             xl(ix,1) = 5d0 * xl2 * xl2
             xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
          else
             call error('pri012','power of L not supported',2)
          end if
       end do

       chi(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
       chi(ipri,2) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*xl(2,0)*xl(3,0)*ex
       chi(ipri,3) = (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(1,0)*xl(3,0)*ex
       chi(ipri,4) = (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*xl(1,0)*xl(2,0)*ex
       chi(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0)&
          +4*al*al*dd(1,iat)**(l(1)+2))*xl(2,0)*xl(3,0)*ex
       chi(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0)&
          +4*al*al*dd(2,iat)**(l(2)+2))*xl(3,0)*xl(1,0)*ex
       chi(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0)&
          +4*al*al*dd(3,iat)**(l(3)+2))*xl(1,0)*xl(2,0)*ex
       chi(ipri,8) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*&
          (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(3,0)*ex
       chi(ipri,9) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*&
          (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*xl(2,0)*ex
       chi(ipri,10)= (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*&
          (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(1,0)*ex

       do ix = 1, 10
          ldopri(ipri,ix) = (abs(chi(ipri,ix))*maxc(ipri) > cutoff_pri)
       enddo
    enddo ! ipri = 1, npri

    ! build the MO avlues at the point
    phi = 0d0
    do ix = 1, 10
       do ipri = 1, m%npri
          if (.not.ldopri(ipri,ix)) cycle
          do imo = 1, m%nmo
             phi(imo,ix) = phi(imo,ix) + m%c(imo,ipri)*chi(ipri,ix)
          enddo
       enddo
    enddo

    ! contribution to the density, etc.
    pr%rho = 0d0
    pr%rhox = 0d0
    pr%rhoxx = 0d0
    pr%bxdm = 0d0
    pr%drho2 = 0d0
    pr%d2rho = 0d0
    pr%taup = 0d0
    gg = 0d0
    hh = 0d0
    if (m%wfntyp == 0) then
       ! closed-shell wavefunction
       do imo = 1, m%nmo
          aocc = m%occ(imo) * 0.5d0
          pr%rho(1) = pr%rho(1) + aocc * phi(imo,1) * phi(imo,1)
          gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
          hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
          pr%taup(1) = pr%taup(1) + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
       enddo
       pr%rho(2) = pr%rho(1)
       pr%rho(0) = pr%rho(1)+pr%rho(2)
       pr%rhox(1,:) = gg
       pr%rhox(2,:) = gg
       pr%rhox(0,:) = pr%rhox(1,:)+pr%rhox(2,:)
       pr%rhoxx(1,:) = hh
       pr%rhoxx(2,:) = hh
       pr%rhoxx(0,:) = pr%rhoxx(1,:)+pr%rhoxx(2,:)
       pr%drho2(1) = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
       pr%drho2(2) = pr%drho2(1)
       pr%drho2(0) = 4 * pr%drho2(1)
       pr%d2rho(1) = hh(1)+hh(2)+hh(3)
       pr%d2rho(2) = pr%d2rho(1)
       pr%d2rho(0) = pr%d2rho(1) + pr%d2rho(2) 
       pr%taup(2) = pr%taup(1)
       pr%taup(0) = pr%taup(1) + pr%taup(2)
       
       ! br hole, one-normalization
       if (pr%rho(1) > small) then
          dsigs = pr%taup(1) - 0.25d0 * pr%drho2(1) / max(pr%rho(1),1d-30)
          quads = (pr%d2rho(1) - 2d0 * dsigs) / 6d0
          call bhole(pr%rho(1),quads,1d0,pr%bxdm(1),pr%alf(1),pr%prefac(1))
          pr%bxdm(2) = pr%bxdm(1)
          pr%alf(2) = pr%alf(1)
          pr%prefac(2) = pr%prefac(1)
          pr%bxdm(0) = pr%bxdm(1)
          pr%alf(0) = pr%alf(1)
          pr%prefac(0) = pr%prefac(1)
       endif
    else if (m%wfntyp == 1) then
       ! open-shell, alpha spin
       nmo1 = (m%nmo + m%mult - 1)/2
       do imo = 1, nmo1
          aocc = m%occ(imo)
          pr%rho(1) = pr%rho(1) + aocc * phi(imo,1) * phi(imo,1)
          gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
          hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
          pr%taup(1) = pr%taup(1) + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
       enddo
       pr%rhox(1,:) = gg
       pr%rhoxx(1,:) = hh

       ! br hole, one-normalization, alpha spin
       if (pr%rho(1) > small) then
          pr%drho2(1) = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
          pr%d2rho(1) = hh(1)+hh(2)+hh(3)
          dsigs = pr%taup(1) - 0.25d0 * pr%drho2(1) / max(pr%rho(1),1d-30)
          quads = (pr%d2rho(1) - 2d0 * dsigs) / 6d0
          call bhole(pr%rho(1),quads,1.d0,pr%bxdm(1),pr%alf(1),pr%prefac(1))
       endif

       ! open-shell, beta spin
       gg = 0d0
       hh = 0d0
       do imo = nmo1+1, m%nmo
          aocc = m%occ(imo)
          pr%rho(2) = pr%rho(2) + aocc * phi(imo,1) * phi(imo,1)
          gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
          hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
          pr%taup(2) = pr%taup(2) + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
       enddo
       pr%rhox(2,:) = gg
       pr%rhoxx(2,:) = hh

       ! br hole, one-normalization, beta spin
       if (pr%rho(2) > small) then
          pr%drho2(2) = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
          pr%d2rho(2) = hh(1)+hh(2)+hh(3)
          dsigs = pr%taup(2) - 0.25d0 * pr%drho2(2) / max(pr%rho(2),1d-30)
          quads = (pr%d2rho(2) - 2d0 * dsigs) / 6d0
          call bhole(pr%rho(2),quads,1d0,pr%bxdm(2),pr%alf(2),pr%prefac(2))
       endif

       ! totals
       pr%rho(0) = pr%rho(1)+pr%rho(2)
       pr%rhox(0,:) = pr%rhox(1,:)+pr%rhox(2,:)
       pr%rhoxx(0,:) = pr%rhoxx(1,:)+pr%rhoxx(2,:)
       pr%drho2(0) = pr%rhox(0,1)*pr%rhox(0,1)+pr%rhox(0,2)*pr%rhox(0,2)+pr%rhox(0,3)*pr%rhox(0,3)
       pr%d2rho(0) = pr%d2rho(1) + pr%d2rho(2) 
       pr%taup(0) = pr%taup(1) + pr%taup(2)

       ! br hole, one-normalizatoin, spin avg
       if (pr%rho(0) > small) then
          dsigs = pr%taup(0) - 0.25d0 * pr%drho2(0) / max(pr%rho(0),1d-30)
          quads = (pr%d2rho(0) - 2d0 * dsigs) / 6d0
          call bhole(0.5d0*pr%rho(0),0.5d0*quads,1d0,pr%bxdm(0),pr%alf(0),pr%prefac(0))
       endif

       ! restricted-open wavefunction

    else if (m%wfntyp == 2) then
       call error("propts_grid","ro not implemented",2)
       ! nmo1 = m%nmo - m%mult + 1
       ! do imo = 1, nmo1
       !    aocc = m%occ(imo) * 0.5d0
       !    pr%rho(2) = pr%rho(2) + aocc * phi(imo,1) * phi(imo,1)
       !    gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
       !    hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
       !    pr%taup(2) = pr%taup(2) + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
       ! enddo
       ! if (pr%rho(2) > small) then
       !    pr%drho2(2) = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
       !    pr%d2rho(2) = hh(1)+hh(2)+hh(3)
       !    dsigs = pr%taup(2) - 0.25d0 * pr%drho2(2) / max(pr%rho(2),1d-30)
       !    quads = (pr%d2rho(2) - 2d0 * dsigs) / 6d0
       !    call bhole(pr%rho(2),quads,1d0,pr%bxdm(2))
       ! endif
       ! pr%rho(1) = pr%rho(2)
       ! do imo = nmo1+1, m%nmo
       !    aocc = m%occ(imo)
       !    pr%rho(1) = pr%rho(1) + aocc * phi(imo,1) * phi(imo,1)
       !    gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
       !    hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
       !    pr%taup(1) = pr%taup(1) + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
       ! enddo
       ! if (pr%rho(1) > small) then
       !    pr%drho2(1) = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
       !    pr%d2rho(1) = hh(1)+hh(2)+hh(3)
       !    dsigs = pr%taup(1) - 0.25d0 * pr%drho2(1) / max(pr%rho(1),1d-30)
       !    quads = (pr%d2rho(1) - 2d0 * dsigs) / 6d0
       !    call bhole(pr%rho(1),quads,1d0,pr%bxdm(1))
       ! endif
    else
       call error("propts_grid","wfn type not implemented",2)
    endif

    if (present(ophi)) then
       ophi = phi(:,1)
    endif

  end subroutine mopoint

  subroutine calc_exdens(m,mesh)
    use meshmod
    use types
    use param

    type(molecule), intent(in) :: m 
    type(tmesh), intent(inout) :: mesh 

    integer :: i, j, ispin, nspin, nmo1, il, iu
    real*8 :: phii(mesh%n), phij(mesh%n), prod(mesh%n), coul(mesh%n)
    real*8 :: fac, exx

    if (.not.didimo.or..not.allocated(mesh%exdens)) &
       call error('calc_exdens','missing phi or matrix alloc',2)

    nspin = 1
    if (m%wfntyp == 1) nspin = 2
    nmo1 = (m%nmo + m%mult - 1)/2

    mesh%exdens = 0d0
    do ispin = 1, nspin
       rewind(iimo)
       if (nspin == 1) then
          il = 1
          iu = m%nmo
       else
          if (ispin == 1) then
             il = 1
             iu = nmo1
          else
             il = nmo1+1
             iu = m%nmo
          end if
       endif
       do i = 1, m%nmo
          read(iimo) phii(1:mesh%n)
          rewind(iimo)
          do j = 1, i
             read(iimo) phij(1:mesh%n)
             if (m%occ(i) < 1d-3 .or. m%occ(j) < 1d-3) cycle
             if (i < il .or. i > iu .or. j < il .or. j > iu) cycle
             prod = phii * phij
             call poiss_grid(m,mesh,0,prod,coul)
             if (nspin == 1) then
                fac = 1d0
             else
                fac = m%occ(i)*m%occ(j)
             end if
             if (i == j) fac = 0.5d0 * fac
             mesh%exdens(:,ispin) = mesh%exdens(:,ispin) -fac * prod * coul
          end do
       end do
       mesh%exdens(:,ispin) = min(mesh%exdens(:,ispin),-3d0*epsilon(1d0))
    end do

    if (nspin == 1) then
       mesh%exdens(:,2) = mesh%exdens(:,1)
       mesh%exdens(:,0) = mesh%exdens(:,1)*2
    else
       mesh%exdens(:,0) = mesh%exdens(:,1)+mesh%exdens(:,2)
    endif

  end subroutine calc_exdens

  !> Calculate the inverse BR hole normalization on the grid
  subroutine calc_xlns(mesh)
    use meshmod
    use types
    use tools
    use param
    implicit none

    type(tmesh), intent(inout) :: mesh 

    integer :: i, j
    real*8 :: quads, dsigs, uxps

    real*8, parameter :: tiny = 1d-20

    if (.not.allocated(mesh%exdens)) &
       call error('calc_xlns','need exdens for xlns',2)
    
    if (.not.allocated(mesh%xlns)) allocate(mesh%xlns(mesh%n,2))
    mesh%xlns = 0d0
    do j = 1, 2
       do i = 1, mesh%n
          dsigs = mesh%tau(i,j) - 0.25d0 * mesh%drho2(i,j) / max(mesh%rho(i,j),1d-30)
          quads = (mesh%d2rho(i,j)-2d0*dsigs) / 6d0
          uxps = -2d0*mesh%exdens(i,j)/max(mesh%rho(i,j),tiny)
          call xlnorm(mesh%rho(i,j),quads,uxps,mesh%xlns(i,j))
       end do
    end do

  end subroutine calc_xlns

end module wfnmod
