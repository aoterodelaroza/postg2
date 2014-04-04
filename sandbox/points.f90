! -*- mode: F90 -*-

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
