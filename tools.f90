! -*- mode: F90 -*-

!> Specific tools that are independent of the mesh or the molecule. 
module tools

  public

  public :: bhole ! calculate the becke-roussel hole parameters
  public :: xlnorm ! inverse BR hole normalization
  public :: xfuncs ! rhs of the BR hole equation and derivative
  
contains
  
  !> Becke-Roussel hole routine
  subroutine bhole(rho,quad,hnorm,b,alf,prefac)
    use param

    real*8, intent(in) :: rho !< Density
    real*8, intent(in) :: quad !< Q-value
    real*8, intent(in) :: hnorm !< Hole normalization
    real*8, intent(out) :: b !< Hole b
    real*8, intent(out) :: alf !< Hole a
    real*8, intent(out) :: prefac !< Hole A

    real*8 :: rhs, x0, shift, x1, x, expo, f, df
    integer :: i

    rhs=third2*(pi*rho/hnorm)**third2*rho/quad
    x0=2.d0
    shift=1.d0
    if(rhs.lt.0.d0)go to 10
    if(rhs.gt.0.d0)go to 20
10  do i=1,16
       x=x0-shift
       call xfuncs(x,rhs,f,df)
       if(f.lt.0.d0)go to 88
       shift=0.1d0*shift
    enddo
    write(iout,1002)
    stop
20  do i=1,16
       x=x0+shift
       call xfuncs(x,rhs,f,df)
       if(f.gt.0.d0)go to 88
       shift=0.1d0*shift
    enddo
    write(iout,1002)
    stop
88  continue
    do i=1,100
       call xfuncs(x,rhs,f,df)
       x1=x-f/df
       if(dabs(x1-x).lt.1.d-10)go to 111
       x=x1
    enddo
    write(iout,1001)
    stop
111 x=x1
    expo=dexp(-x)
    prefac=rho/expo
    alf=(8.d0*pi*prefac/hnorm)**third
    b=x/alf
    return
1001 format(' ','bhole: newton algorithm fails to converge!')
1002 format(' ','bhole: newton algorithm fails to initialize!')
  end subroutine bhole

  !> Calculate the inverse BR hole normalization. From numol.
  subroutine xlnorm(rho,quad,uxpos,xlnrm)
    use param
    implicit none

    real*8, intent(in) :: rho, quad, uxpos
    real*8, intent(out) :: xlnrm

    real*8 :: rhs, x0, shift, x, f, df, x1, alf, a
    integer :: i
    logical :: found

    if (rho < 1.d-10) then
       xlnrm=1.D0
       return
    end if
    rhs=4.d0*pi/3.d0*rho*rho/quad/uxpos
    x0=2.D0
    shift=1.D0
    found = .false.
    if (rhs < 0.d0) then
       do i = 1, 16
          x = x0 - shift
          call xlfuns(x,rhs,f,df)
          if (f < 0.D0) then
             found = .true.
             exit
          end if
          shift=0.1d0*shift
       end do
       if (.not.found) &
          call error('xlnorm','newton algorithm failed to initialize',2)
    else
       do i=1, 16
          x=x0+shift
          call xlfuns(x,rhs,f,df)
          if(f > 0.D0) then
             found = .true.
             exit
          endif
          shift=0.1d0*shift
       end DO
       if (.not.found) &
          call error('xlnorm','newton algorithm failed to initialize',2)
    endif
    found = .false.
    do i = 1, 100
       call xlfuns(x,rhs,f,df)
       x1=x-f/df
       if(dabs(x1-x) < 1.d-10) then
          found = .true.
          exit
       end if
       x=x1
    end do
    if (.not.found) &
       call error('xlnorm','newton algorithm failed to converge',2)
    x=x1
    alf=dsqrt(6.d0*quad*x/rho/(x-2.d0))
    a=rho*exp(x)
    xlnrm=min(8.d0*pi*a/alf**3,2.d0)

  end subroutine xlnorm

  !< RHS of the hole equation, and derivative. 
  subroutine xfuncs(x,rhs,f,df)
    real*8, intent(in) :: x !< x-value
    real*8, intent(in) :: rhs !< Right-hand side
    real*8, intent(out) :: f !< Value of the function
    real*8, intent(out) :: df !< Value of the function derivative

    real*8 :: expo23

    expo23=dexp(-2.d0/3.d0*x)
    f = x*expo23/(x-2.d0) - rhs
    df=2.d0/3.d0*(2.d0*x-x**2-3.d0)/(x-2.d0)**2*expo23
  end subroutine xfuncs

  !< RHS of the BR hole equation and derivative
  subroutine xlfuns(x,rhs,f,df)
    implicit none

    real*8, intent(in) :: x, rhs
    real*8, intent(out) :: f, df

    real*8 :: expo, bot

    expo=exp(x)
    bot=(x-2.d0)*(expo-1.d0-0.5d0*x)
    f = x*x/bot - rhs
    df=4.d0*x-(4.d0*x-3.d0*x*x+x**3)*expo
    df=df/bot**2
  end subroutine xlfuns

end module tools
