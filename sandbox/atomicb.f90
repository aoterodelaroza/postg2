! -*- mode: F90 -*-

!> Calculate the atomic average of BR hole's b
subroutine sandbox_atomicb(mol)
  use wfnmod
  use meshmod
  use tools_math
  use types
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule

  integer :: nr, nang 

  real*8 :: rads(1400), wrads(1400), x(3), alf, prefac
  real*8 :: rq(1400), rqq(1400), d2r(1400,7), trid2r(1400,3)
  real*8 :: xang(302), yang(302), zang(302), wang(302)
  real*8 :: rmid, bb1, bb2, rho, taup, drho2, d2rho, dsigs, quads
  real*8 :: alf1
  type(props) :: pr
  integer :: i, ir
  logical :: mask(mprops)

  ! radial mesh
  if (mol%z(1) <= 2) then
     nr = 200
  elseif (mol%z(1) <= 10) then
     nr = 400
  elseif (mol%z(1) <= 18) then
     nr = 600
  elseif (mol%z(1) <= 36) then
     nr = 800
  elseif (mol%z(1) <= 54) then
     nr = 1000
  elseif (mol%z(1) <= 86) then
     nr = 1200
  elseif (mol%z(1) <= 94) then
     nr = 1400
  else
     call error('atomin','atomic number out of range',2)
  endif
  rmid = 1d0 / real(mol%z(1),8)**third
  call rmesh(nr,rmid,rads,rq,rqq,wrads,d2r,trid2r)

  ! angular mesh
  nang = 302
  call ld0302(xang,yang,zang,wang,nang)
  
  mask = .false.
  mask(1:19) = .true.
  write (iout,'("# r bavg1 rho")')
  aloop: do ir = 1, nr
     bb1 = 0d0
     alf1 = 0d0
     rho = 0d0
     taup = 0d0
     drho2 = 0d0
     d2rho = 0d0
     do i = 1, nang
        x = (/xang(i),yang(i),zang(i)/) * rads(ir)
        call propts_point(mol,x,pr,mask)

        ! first average
        ! bb1 = bb1 + wang(i) * pr%bxdm(0)
        bb1 = bb1 + wang(i) * (pr%rho(1)*pr%bxdm(1)+pr%rho(2)*pr%bxdm(2)) / max(pr%rho(0),1d-30)
        alf1 = alf1 + wang(i) * (pr%rho(1)*pr%alf(1)+pr%rho(2)*pr%alf(2)) / max(pr%rho(0),1d-30)

        ! for second average
        rho = rho + wang(i) * 0.5d0 * pr%rho(0)
        taup = taup + wang(i) * 0.5d0 * pr%taup(0)
        drho2 = drho2 + wang(i) * 0.25d0 * pr%drho2(0)
        d2rho = d2rho + wang(i) * 0.5d0 * pr%d2rho(0)
     end do
     if (rho > 1d-5) then
        dsigs = taup - 0.25d0 * drho2 / max(rho,1d-30)
        quads = (d2rho - 2d0 * dsigs) / 6d0
        call bhole(rho,quads,1d0,bb2,alf,prefac)
        write (iout,'(1p,5(E20.12,X))') rads(ir), bb1, alf1
     else
        exit
     endif
  end do aloop

  ! spinpol
!  write (iout,'("# r bup bdn rhoup rhodn")')
!  do ir = 1, nr
!     bb = 0d0
!     rr = 0d0
!     do i = 1, nang
!        x = (/xang(i),yang(i),zang(i)/) * rads(ir)
!        call propts_point(mol,x,pr)
!        bb = bb + wang(i) * pr%bxdm
!        rr = rr + wang(i) * pr%rho
!     end do
!     if (rr(1)+rr(2) > 1d-6) then
!        write (iout,'(1p,5(E20.12,X))') rads(ir), bb, rr
!     endif
!  end do

101 continue

  return
999 call error("sandbox_points","postg2 file.wfx points iprop << x0 y0 z0 ",2)

end subroutine sandbox_atomicb
