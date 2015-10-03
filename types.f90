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

!> Module for data types and memory management
module types
  use param, only:mline
  implicit none
  public

  public :: propts_select ! extract property number n from a property object
  public :: mask_all ! return an all-true properties mask
  public :: mask_none ! return an all-false properties mask
  public :: mesh_allocate ! allocate the mesh arrays using a mask

  integer, parameter :: mprops = 22
  !> Molecular mesh
  type tmesh
     integer :: n
     integer, allocatable :: ktab(:) !< mesh indices for nuclei
     real*8, allocatable :: w(:) !< mesh weights
     real*8, allocatable :: wnuc(:) !< mesh weights
     real*8, allocatable :: x(:,:) !< mesh coordinates
     real*8, allocatable :: promol(:) !< promolecular density on the grid
     real*8, allocatable :: bpromol(:) !< b on the grid
     real*8, allocatable :: rho(:,:) !< density (up and down) 
     real*8, allocatable :: drho2(:,:) !< norm of the gradient 
     real*8, allocatable :: d2rho(:,:) !< laplacian 
     real*8, allocatable :: tau(:,:) !< kinetic energy density
     real*8, allocatable :: bxdm(:,:) !< b BR hole parameter, for XDM (1-norm)
     real*8, allocatable :: alf(:,:)  !< a BR hole parameter (1-norm)
     real*8, allocatable :: prefac(:,:) !< A BR hole parameter (1-norm)
     real*8, allocatable :: vel(:) !< electrostatic potential
     real*8, allocatable :: exdens(:,:) !< exchange energy density
     real*8, allocatable :: xlns(:,:) !< inverse BR hole normalization
     logical :: isthere(mprops) = .false. !< true if a property is calculated
  end type tmesh
  !> Molecule
  type molecule
     character(mline) :: name !< name of the molecule
     integer :: n !< number of atoms
     real*8, allocatable :: x(:,:) !< atomic coordinates
     integer, allocatable :: z(:) !< atomic numbers
     integer :: nmo !< number of molecular orbitals
     integer :: npri !< number of primitives
     integer, allocatable :: icenter(:) !< labels for primitive centers
     integer, allocatable :: itype(:) !< labels for primitive types
     real*8, allocatable :: e(:) !< exponents
     real*8, allocatable :: eps(:) !< molecular orbital energies
     real*8, allocatable :: occ(:) !< occupation numbers (orbital)
     real*8, allocatable :: c(:,:) !< orbital coefficients
     integer :: mult !< multiplicity
     integer :: wfntyp  !< 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
     real*8 :: nelec !< number of electrons
     real*8 :: charge !< charge
     real*8 :: escf !< energy
     logical :: useecp !< are we using ECPs?
  end type molecule
  !> Properties
  type props
     real*8 :: rho(0:2)     !< density (total, up, dn)
     real*8 :: rhox(0:2,3)  !< density gradient (total, up, dn)
     real*8 :: rhoxx(0:2,3) !< density hessian diagonal (total, up, dn)
     real*8 :: drho2(0:2)   !< square of gradient (total, up, dn)
     real*8 :: d2rho(0:2)   !< laplacian (total, up, dn)
     real*8 :: taup(0:2)    !< kinetic energy density (total, up, dn)
     real*8 :: bxdm(0:2)    !< BR hole b, with one-normalization (spin-avg, up, dn)
     real*8 :: alf(0:2)     !< BR hole a, with one-normalization (spin-avg, up, dn)
     real*8 :: prefac(0:2)  !< BR hole A, with one-normalization (spin-avg, up, dn)
     real*8 :: vel          !< electrostatic potential
  end type props

contains

  !> Extract property number n from properties list n
  function propts_select(pr,n) result(x)
    use param, only: error

    type(props), intent(in) :: pr !< Properties list
    integer, intent(in) :: n !< Property to extract
    real*8 :: x !< Property value

    if (n == 1) then
       x = pr%rho(0)
    elseif (n == 2) then
       x = pr%rho(1)
    elseif (n == 3) then
       x = pr%rho(2)
    elseif (n == 4) then
       x = pr%rho(1)-pr%rho(2)
    elseif (n == 5) then
       x = pr%drho2(0)
    elseif (n == 6) then
       x = pr%drho2(1)
    elseif (n == 7) then
       x = pr%drho2(2)
    elseif (n == 8) then
       x = pr%drho2(1)-pr%drho2(2)
    elseif (n == 9) then
       x = pr%d2rho(0)
    elseif (n == 10) then
       x = pr%d2rho(1)
    elseif (n == 11) then
       x = pr%d2rho(2)
    elseif (n == 12) then
       x = pr%d2rho(1)-pr%d2rho(2)
    elseif (n == 13) then
       x = pr%taup(0)
    elseif (n == 14) then
       x = pr%taup(1)
    elseif (n == 15) then
       x = pr%taup(2)
    elseif (n == 16) then
       x = pr%taup(1)-pr%taup(2)
    elseif (n == 17) then
       x = pr%bxdm(0)
    elseif (n == 18) then
       x = pr%bxdm(1)
    elseif (n == 19) then
       x = pr%bxdm(2)
    elseif (n == 20) then
       x = pr%vel
    else
       call error('propts_select','property not found',2)
    end if

  end function propts_select

  !> Return a properties mask with all properties active
  function mask_all()
    logical :: mask_all(mprops)
    mask_all = .true.
  end function mask_all

  !> Return a properties mask with no properties active
  function mask_none()
    logical :: mask_none(mprops)
    mask_none = .false.
  end function mask_none

  !> Allocate the arrays in the mesh using a mask
  subroutine mesh_allocate(mesh,mask)
    use param, only: error
    type(tmesh), intent(inout) :: mesh !< Mesh to be allocated
    logical, intent(in) :: mask(mprops) !< Input mask

    integer :: istat

    istat = 0
    ! density
    if (any(mask(1:4))) then
       if (.not.allocated(mesh%rho)) allocate(mesh%rho(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for rho',2)
    end if
    ! gradient
    if (any(mask(5:8))) then
       if (.not.allocated(mesh%drho2)) allocate(mesh%drho2(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for grad',2)
    end if
    ! laplacian
    if (any(mask(9:12))) then
       if (.not.allocated(mesh%d2rho)) allocate(mesh%d2rho(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for lap',2)
    end if
    ! ked
    if (any(mask(13:16))) then
       if (.not.allocated(mesh%tau)) allocate(mesh%tau(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for tau',2)
    end if
    ! BR b
    if (mask(17)) then
       if (.not.allocated(mesh%bxdm)) allocate(mesh%bxdm(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for bxdm',2)
    end if
    ! BR a
    if (mask(18)) then
       if (.not.allocated(mesh%alf)) allocate(mesh%alf(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for alf',2)
    end if
    ! BR A
    if (mask(19)) then
       if (.not.allocated(mesh%prefac)) allocate(mesh%prefac(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for prefac',2)
    end if
    ! Electrostatic potential
    if (mask(20)) then
       if (.not.allocated(mesh%vel)) allocate(mesh%vel(mesh%n),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for vel',2)
    end if
    ! Exchange potential
    if (mask(21)) then
       if (.not.allocated(mesh%exdens)) allocate(mesh%exdens(mesh%n,0:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for exdens',2)
    end if
    ! Inverse BR hole normalization
    if (mask(22)) then
       if (.not.allocated(mesh%xlns)) allocate(mesh%xlns(mesh%n,1:2),stat=istat)
       if (istat /= 0) call error('propts_grid','could not allocate memory for xlns',2)
    end if

  end subroutine mesh_allocate

end module types
