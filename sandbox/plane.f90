! -*- mode: F90 -*-

!> Calculate properties on a plane
subroutine sandbox_plane(mol,mesh)
  use types
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule
  type(tmesh), intent(inout) :: mesh !< Molecular mesh

  character*(mline) :: line
  real*8 :: x0(3), x1(3), x2(3), x(3), pro, d1, d2
  integer :: iprop, n1, n2, i, j
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
  read (line,*,end=999,err=999) x2(1)
  call getarg(11,line)
  read (line,*,end=999,err=999) x2(2)
  call getarg(12,line)
  read (line,*,end=999,err=999) x2(3)
  call getarg(13,line)
  read (line,*,end=999,err=999) n1
  call getarg(14,line)
  read (line,*,end=999,err=999) n2

  mask = .false.
  mask(iprop) = .true.

  write (iout,'("# x y z u v value")') 
  do i = 0, n1-1
     do j = 0, n2-1
        x = x0 + (x1-x0) * dble(i) / (n1-1) + (x2-x0) * dble(j) / (n2-1)
        d1 = sqrt(dot_product(x1-x0,x1-x0)) * dble(i)/(n1-1)
        d2 = sqrt(dot_product(x2-x0,x2-x0)) * dble(j)/(n2-1)
        call propts_point(mol,x,pr,mask,mesh)
        pro  = propts_select(pr,iprop)
        write (iout,'(5(F16.10,X),1p,E20.12)') x, d1, d2, pro
     end do
     write (iout,*)
  end do

101 continue

  return
999 call error("sandbox_plane","postg2 file.wfx plane iprop x0 y0 z0 x1 y1 z1 x2 y2 z2 n1 n2",2)

end subroutine sandbox_plane
