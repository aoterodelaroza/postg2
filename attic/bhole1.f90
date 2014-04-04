  !< new hole model, with b given by the promolecular average
  subroutine bhole1(rho,b,alf,prefac)
    use param

    real*8, intent(in) :: rho
    real*8, intent(in) :: b
    real*8, intent(out) :: alf, prefac

    integer, parameter :: nmax = 100
    real*8 :: fa, fpa, alf1
    real*8 :: alfa, alfb, a, fb, fx

    integer :: i

    ! bisection
    if (b < 1.5d0) then
       alfa = 0d0
       alfb = 3d0 / b
    else
       alfa = 3d0 / b 
       alfb = 10d0
    endif
    fa = alfa**3 * exp(-b*alfa) - 8d0 * pi * rho
    fb = alfb**3 * exp(-b*alfb) - 8d0 * pi * rho
    do while(abs(alfa-alfb) > 1d-10)
       alf = 0.5d0 * (alfa + alfb)
       fx = alf**3 * exp(-b*alf) - 8d0 * pi * rho
       if (fa*fx > 0d0) then
          alfa = alf
          fa = fx
       else
          alfb = alf
          fb = fx
       endif
    end do

    ! ! newton
    ! alf1 = 0d0
    ! alf = 3d0 / 2d0 / b
    ! do i = 1, nmax
    !    fa = alf**3 * exp(-b*alf) - 8d0 * pi * rho
    !    fpa = (3-b*alf) * alf**2 * exp(-b*alf)
    !    if (abs(fpa) < 1d-5) call error("bhole1","fpa = 0",2)
    !    alf1 = alf
    !    alf = alf - fa / fpa
    !    if (abs(alf-alf1) < 1d-10) exit
    !    if (alf < 0 .or. alf > 3d0/b) call error("bhole1","alf out of bounds",2)
    ! end do
    ! if (i == nmax) call error("bhole1","newton failed to converge",2)

    prefac = alf**3 / 8d0 / pi

  end subroutine bhole1


