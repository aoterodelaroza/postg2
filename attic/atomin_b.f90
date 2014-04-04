  subroutine atomin_b(m,mesh)
    use tools_math, only: spline
    use promolmod
    implicit none
    ! atomic b averages

    type(molecule), intent(in) :: m
    type(tmesh), intent(inout) :: mesh

    character*(mline) :: postg_home, afile
    real*8 :: rmid, rdum(8), h, q, x(3), r, dq, arho, ab, rdata
    integer :: intq, ilast
    real*8 :: blast, rlast, bend, r0, r1
    real*8, allocatable, dimension(:) :: a2, b2, c2, f2
    integer :: i, j, l, ndata, idum1, idum2, is, istat, kk, nn
    character*2 :: compar
    logical :: ok
    integer :: isenv

    if (.not.allocated(mesh%bpromol)) allocate(mesh%bpromol(mesh%n),stat=istat)
    if (istat /= 0) call error('atominb','could not allocate bpromol',2)
    mesh%bpromol = 0d0

    do i = 1, m%n
       if (m%z(i) < 1) cycle

       ! number of data points
       rmid = 1d0 / m%z(i)**third
       if (m%z(i) <= 2) then
          ndata = 200
       elseif (m%z(i) <= 10) then
          ndata = 400
       elseif (m%z(i) <= 18) then
          ndata = 600
       elseif (m%z(i) <= 36) then
          ndata = 800
       elseif (m%z(i) <= 54) then
          ndata = 1000
       elseif (m%z(i) <= 86) then
          ndata = 1200
       elseif (m%z(i) <= 94) then
          ndata = 1400
       else
          call error('atomin','atomic number out of range',2)
       endif
       h = 1d0 / (ndata+1)

       ! complete the b array by linear extrapolation
       if (abs(bavg(1,m%z(i))) < 1d-10) call error('atominb','atom not found',2)
       do j = 1, ndata
          if (abs(bavg(j,m%z(i))) > 1d-10) then
             q = h * j
             ilast = j
             rlast = rmid * q / (1d0 - q)
             blast = bavg(j,m%z(i))
          end if
       end do

       ! spline for rho
       allocate(f2(0:ndata),a2(0:ndata),b2(0:ndata),c2(0:ndata))
       if (istat /= 0) call error('atomin','could not allocate memory for f2, a2, b2, c2',2)
       f2 = 0d0
       f2(1:ndata) = ftot(1:ndata,m%z(i))
       call spline(h,f2,a2,b2,c2,ndata,0.d0)

       ! interpolate on the mesh
       do kk = 1, mesh%n
          x = mesh%x(:,kk) - m%x(:,i)
          r = sqrt(x(1)**2+x(2)**2+x(3)**2)

          ! spline interpolation for rho
          q = r / (r + rmid)
          intq = int((ndata+1) * q)
          dq = q - intq * h
          arho = abs((f2(intq)+dq*(a2(intq)+dq*(b2(intq)+dq*c2(intq)))))/r**2

          ! linear interpolation for b
          r0 = intq*h * rmid / (1 - intq*h)
          r1 = (intq+1)*h * rmid / (1 - (intq+1)*h)
          if (intq < ilast) then
             ab = bavg(intq,m%z(i)) + (r-r0) * (bavg(intq+1,m%z(i))-bavg(intq,m%z(i))) / (r1-r0)
          else
             ab = blast + (r - rlast)
          endif

          ! xxxx
          ! mesh%bpromol(kk) = mesh%bpromol(kk) + ab * arho
          mesh%bpromol(kk) = mesh%bpromol(kk) + r * arho
       enddo
       deallocate(f2,a2,b2,c2)
    end do
    mesh%bpromol = mesh%bpromol / mesh%promol

  end subroutine atomin_b

