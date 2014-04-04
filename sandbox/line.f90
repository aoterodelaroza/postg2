! -*- mode: F90 -*-

!> Calculate properties on a line
subroutine sandbox_line(mol,mesh)
  use types
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule 
  type(tmesh), intent(inout) :: mesh !< Molecular mesh

  character*(mline) :: line
  real*8 :: x0(3), x1(3), x(3), pro
  integer :: iprop, npts, i
  type(props) :: pr
  logical :: mask(mprops)

  call getarg(3,line)
  read (line,*,end=999,err=999) iprop
  call getarg(4,line)
  read (line,*,end=999,err=999) x0(1)
  call getarg(5,line)
  read (line,*,end=999,err=999) x0(2)
  call getarg(6,line)
  read (line,*,end=999,err=999) x0(3)
  call getarg(7,line)
  read (line,*,end=999,err=999) x1(1)
  call getarg(8,line)
  read (line,*,end=999,err=999) x1(2)
  call getarg(9,line)
  read (line,*,end=999,err=999) x1(3)
  call getarg(10,line)
  read (line,*,end=999,err=999) npts

  mask = .false.
  mask(iprop) = .true.

  do i = 0, npts-1
     x = x0 + (x1-x0) * dble(i) / (npts-1)
     call propts_point(mol,x,pr,mask,mesh)
     pro  = propts_select(pr,iprop)
     write (iout,'(3(F16.10,X),1p,E20.12)') x, pro
  end do

101 continue

  return
999 call error("sandbox_line","postg2 file.wfx line iprop x0 y0 z0 x1 y1 z1 npts",2)

end subroutine sandbox_line
