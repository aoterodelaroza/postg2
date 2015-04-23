! -*- mode: F90 -*-

!> Interface to the sandbox routines
module sandbox
  use wfnmod
  use io
  use types
  use tools
  use param
  implicit none

  public
  
  public :: driver_sandbox ! parse the commands, initialize, and call the sandbox routine

contains

  !> Parse the command line and call the appropriate sandbox routine.
  subroutine driver_sandbox(line,mol,wfnext,mesh)

    character*(mline), intent(inout) :: line !< Input line
    type(molecule), intent(inout) :: mol !< Molecule
    character*11, intent(inout) :: wfnext !< Wavefunction extension
    type(tmesh), intent(inout) :: mesh !< Molecular mesh
    logical :: mask(mprops)
    
    ! read the first command
    call getarg(2,line)
    if (trim(adjustl(lower(line))) == "info") then
       call popwfninfo(mol,wfnext)
    elseif (trim(adjustl(lower(line))) == "xdm") then
       mask = .false.
       mask(1:17) = .true.
       call popwfninfo(mol,wfnext)
       call propts_grid(mol,mesh,mask)
       call sandbox_xdm(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "alpha") then
       mask = .false.
       mask(1:17) = .true.
       call propts_grid(mol,mesh,mask)
       call sandbox_alpha(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "xdmc9") then
       mask = .false.
       mask(1:17) = .true.
       call propts_grid(mol,mesh,mask)
       call sandbox_xdmc9(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "xyz") then
       call sandbox_xyz(mol)
    elseif (trim(adjustl(lower(line))) == "points" .or. trim(adjustl(lower(line))) == "point") then
       call sandbox_points(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "line") then
       call sandbox_line(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "plane") then
       call sandbox_plane(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "cube") then
       call sandbox_cube(mol,mesh)
    elseif (trim(adjustl(lower(line))) == "atomicb") then
       call sandbox_atomicb(mol)
#ifdef HAVE_LIBXC
    elseif (trim(adjustl(lower(line))) == "cube_libxc" .or.&
       trim(adjustl(lower(line))) == "libxc_cube") then
       call sandbox_cubelibxc(mol,mesh)
#endif
    elseif (trim(adjustl(lower(line))) == "energy") then
       mask = .false.
       mask(1:19) = .true.
       call propts_grid(mol,mesh,mask)
       call sandbox_energy(mol,mesh)
    else
       write (iout,*) "command not found!"
    endif
  end subroutine driver_sandbox

#include "sandbox/xdm.f90"
#include "sandbox/xdmc9.f90"
#include "sandbox/xyz.f90"
#include "sandbox/atomicb.f90"
#include "sandbox/points.f90"
#include "sandbox/line.f90"
#include "sandbox/plane.f90"
#include "sandbox/cube.f90"
#include "sandbox/cube_libxc.f90"
#include "sandbox/energy.f90"
#include "sandbox/alpha.f90"

end module sandbox
