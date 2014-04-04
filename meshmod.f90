! -*- mode: F90 -*-

!> Manipulate the molecular integration mesh
module meshmod
  implicit none

  private
  public :: genmesh ! generate the mesh from mol (points, nuc-w, int-w)
  public :: rmesh ! generate a radial integration mesh
  public :: poiss_grid ! solve poisson's equation on the grid

contains

  !> Generate the mesh for the input molecule. Mesh points, nuclear weights and 
  !> integration weights.
  function genmesh(m) result(mesh)
    use tools_math
    use types
    use param
    implicit none

    type(molecule), intent(in) :: m !< Input molecule
    type(tmesh) :: mesh !< Output mesh

    real*8 :: rr(m%n,m%n), rmid, r, r1, r2, hypr, vp0, vpsum, vpi
    integer :: i, j, k, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    real*8, allocatable :: rq(:), rqq(:), d2r(:,:), trid2r(:,:)
    integer :: nr, nang, nangj, nrj, ir, il, istat
    real*8 :: cutoff(m%n,m%n), x(3)

    ! interatomic distances
    rr = 0d0
    do i = 1, m%n
       do j = i+1, m%n
          rr(i,j) = sqrt((m%x(1,i)-m%x(1,j))**2+(m%x(2,i)-m%x(2,j))**2+(m%x(3,i)-m%x(3,j))**2)
          rr(j,i) = rr(i,j)
       enddo
    enddo

    ! allocate space for the mesh
    mesh%n = 0
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       mesh%n = mesh%n + z2nr(m%z(i)) * z2nang(m%z(i))
    enddo
    allocate(mesh%w(mesh%n),mesh%wnuc(mesh%n),mesh%x(3,mesh%n),mesh%ktab(m%n))

    ! no parallelization to have the grid nodes in order
    kk = 0
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       mesh%ktab(i) = kk
       ! radial mesh
       nr = z2nr(m%z(i))
       nang = z2nang(m%z(i))
       allocate(rads(nr),wrads(nr),stat=istat)
       allocate(rq(nr),rqq(nr),d2r(nr,7),trid2r(nr,3))
       if (istat /= 0) call error('genmesh','could not allocate memory for radial meshes',2)
       rmid = 1d0/real(m%z(i),8)**third
       call rmesh(nr,rmid,rads,rq,rqq,wrads,d2r,trid2r)

       ! angular mesh
       nang = z2nang(m%z(i))
       allocate(xang(nang),yang(nang),zang(nang),wang(nang),stat=istat)
       if (istat /= 0) call error('readwfn','could not allocate memory for angular meshes',2)
       if (nang == 6) then
          call ld0006(xang,yang,zang,wang,nang)
       elseif (nang == 14) then
          call ld0014(xang,yang,zang,wang,nang)
       elseif (nang == 26) then
          call ld0026(xang,yang,zang,wang,nang)
       elseif (nang == 38) then
          call ld0038(xang,yang,zang,wang,nang)
       elseif (nang == 50) then
          call ld0050(xang,yang,zang,wang,nang)
       elseif (nang == 74) then
          call ld0074(xang,yang,zang,wang,nang)
       elseif (nang == 86) then
          call ld0086(xang,yang,zang,wang,nang)
       elseif (nang == 110) then
          call ld0110(xang,yang,zang,wang,nang)
       elseif (nang == 146) then
          call ld0146(xang,yang,zang,wang,nang)
       elseif (nang == 170) then
          call ld0170(xang,yang,zang,wang,nang)
       elseif (nang == 194) then
          call ld0194(xang,yang,zang,wang,nang)
       elseif (nang == 230) then
          call ld0230(xang,yang,zang,wang,nang)
       elseif (nang == 266) then
          call ld0266(xang,yang,zang,wang,nang)
       elseif (nang == 302) then
          call ld0302(xang,yang,zang,wang,nang)
       elseif (nang == 350) then
          call ld0350(xang,yang,zang,wang,nang)
       elseif (nang == 434) then
          call ld0434(xang,yang,zang,wang,nang)
       elseif (nang == 590) then
          call ld0590(xang,yang,zang,wang,nang)
       elseif (nang == 770) then
          call ld0770(xang,yang,zang,wang,nang)
       elseif (nang == 974) then
          call ld0974(xang,yang,zang,wang,nang)
       elseif (nang == 1202) then
          call ld1202(xang,yang,zang,wang,nang)
       elseif (nang == 1454) then
          call ld1454(xang,yang,zang,wang,nang)
       elseif (nang == 1730) then
          call ld1730(xang,yang,zang,wang,nang)
       elseif (nang == 2030) then
          call ld2030(xang,yang,zang,wang,nang)
       elseif (nang == 2354) then
          call ld2354(xang,yang,zang,wang,nang)
       elseif (nang == 2702) then
          call ld2702(xang,yang,zang,wang,nang)
       elseif (nang == 3074) then
          call ld3074(xang,yang,zang,wang,nang)
       elseif (nang == 3470) then
          call ld3470(xang,yang,zang,wang,nang)
       elseif (nang == 3890) then
          call ld3890(xang,yang,zang,wang,nang)
       elseif (nang == 4334) then
          call ld4334(xang,yang,zang,wang,nang)
       elseif (nang == 4802) then
          call ld4802(xang,yang,zang,wang,nang)
       elseif (nang == 5294) then
          call ld5294(xang,yang,zang,wang,nang)
       elseif (nang == 5810) then
          call ld5810(xang,yang,zang,wang,nang)
       else
          call error("genmesh","unknown value of nang",2)
       end if

       ! 3d mesh, parallelize over radial shells
       do ir = 1, nr
          r = rads(ir)
          do il = 1, nang
             x = m%x(:,i) + r * (/xang(il),yang(il),zang(il)/)
             do j = 2, m%n
                if (m%z(j) < 1) cycle
                do k = 1, j-1
                   if (m%z(k) < 1) cycle
                   r1 = sqrt((x(1)-m%x(1,j))**2+(x(2)-m%x(2,j))**2+(x(3)-m%x(3,j))**2)
                   r2 = sqrt((x(1)-m%x(1,k))**2+(x(2)-m%x(2,k))**2+(x(3)-m%x(3,k))**2)
                   hypr = (r1-r2) / rr(j,k)
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   cutoff(j,k) = (1d0-hypr) / 2d0
                   cutoff(k,j) = (1d0+hypr) / 2d0
                enddo
                cutoff(j,j) = 1d0
             enddo
             cutoff(1,1) = 1d0
             vp0 = 1d0
             vpsum = 0d0
             do j = 1, m%n
                if (m%z(j) < 1) cycle
                vp0=vp0*cutoff(i,j)
                vpi=1d0
                do k = 1, m%n
                   if (m%z(k) < 1) cycle
                   vpi = vpi * cutoff(j,k)
                enddo
                vpsum = vpsum + vpi
             enddo
             kk = kk + 1
             mesh%wnuc(kk) = vp0/vpsum
             mesh%w(kk) = vp0/vpsum * wrads(ir) * wang(il)
             mesh%x(:,kk) = x
          enddo
       enddo
       deallocate(rads,wrads,xang,yang,zang,wang,rq,rqq,d2r,trid2r)
    enddo

    ! initialize the flags
    mesh%isthere = .false.

  end function genmesh

  !> Calculate the radial integration mesh, weights, and derivatives
  !> of the variable r with respect to q
  !> The q-mesh is uniform on the interval (0,+1). The transformation is 
  !> r = rmid * q / (1-q). Also, 3- and 7- point finite difference 
  !> matrices for d2(wrt)r on q-grid.
  subroutine rmesh(n,rmid,r,rq,rqq,wintr,d2r,trid2r)
    use param
    implicit none

    integer, intent(in) :: n !< Number of points
    real*8, intent(in) :: rmid !< Mid-point
    real*8, intent(out) :: r(n) !< Radial points
    real*8, intent(out) :: rq(n) 
    real*8, intent(out) :: rqq(n) 
    real*8, intent(out) :: wintr(n) !< Radial weights
    real*8, intent(out) :: d2r(n,7)
    real*8, intent(out) :: trid2r(n,3)

    real*8 :: h, q, c1, c2, cr1, cr2, a1, a2
    integer :: i, j, k, ii
    real*8 :: fd1(-3:3,-3:3), fd2(-3:3,-3:3)

    data ((fd1(I,J),J=-3,3),I=-3,3)/&
          0.D0,   0.D0,   -6.D0,   -20.D0,  36.D0,  -12.D0,  2.D0,&
          0.D0,   6.D0,  -60.D0,   -40.D0, 120.D0,  -30.D0,  4.D0,&
        -12.D0, 108.D0, -540.D0,     0.D0, 540.D0, -108.D0, 12.D0,&
        -12.D0, 108.D0, -540.D0,     0.D0, 540.D0, -108.D0, 12.D0,&
        -12.D0, 108.D0, -540.D0,     0.D0, 540.D0, -108.D0, 12.D0,&
         -4.D0,  30.D0, -120.D0,    40.D0,  60.D0,   -6.D0,  0.D0,&
         -2.D0,  12.D0,  -36.D0,    20.D0,   6.D0,    0.D0,  0.D0/
    data ((FD2(I,J),J=-3,3),I=-3,3)/&
          0.D0,   0.D0,   11.D0,   -20.D0,   6.D0,    4.D0, -1.D0,&
          0.D0,  -5.D0,   80.D0,  -150.D0,  80.D0,   -5.D0,  0.D0,&
          4.D0, -54.D0,  540.D0,  -980.D0, 540.D0,  -54.D0,  4.D0,&
          4.D0, -54.D0,  540.D0,  -980.D0, 540.D0,  -54.D0,  4.D0,&
          4.D0, -54.D0,  540.D0,  -980.D0, 540.D0,  -54.D0,  4.D0,&
          0.D0,  -5.D0,   80.D0,  -150.D0,  80.D0,   -5.D0,  0.D0,&
         -1.D0,   4.D0,    6.D0,   -20.D0,  11.D0,    0.D0,  0.D0/

    h = 1d0/real(n+1,8)
    do i=1,n
       q=h*i
       r(i) = rmid * q / (1.d0-q)
       rq(i)=rmid/(1d0-q)**2
       rqq(i)=2d0*rmid/(1.d0-q)**3
       wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
    enddo

    ! construct second derivative matrix
    do i=1,n
       if(i==1.or.i==n)then
          c1=1.d0/24.d0/h
          c2=1.d0/12.d0/h**2
       else if(i==2.or.i==n-1)then
          c1=1.d0/120.d0/h
          c2=1.d0/ 60.d0/h**2
       else
          c1=1.d0/720.d0/h
          c2=1.d0/360.d0/h**2
       end if
       cr2=c2/rq(i)**2
       cr1=-c1*rqq(i)/rq(i)**3
       ii=0
       if(i.le.3)ii=i-4
       if(i.ge.n-2)ii=i-n+3
       do k=1,7
          d2r(i,k)=cr1*fd1(ii,k-4)+cr2*fd2(ii,k-4)
       end do
    end do

    do i=1,n
       a2=1.d0/rq(i)**2/h**2
       a1=-0.5d0*rqq(i)/rq(i)**3/h
       trid2r(i,1)=a2-a1
       trid2r(i,2)=-2.d0*a2
       trid2r(i,3)=a2+a1
    end do

  end subroutine rmesh

  !> Solve Poisson's equation on the grid. Adapted from numol.
  subroutine poiss_grid(mol,mesh,itype,rho0,vel0)
    use types
    use param
    use tools_math
    use tools_linpack

    type(tmesh), intent(inout) :: mesh
    type(molecule), intent(in) :: mol
    integer, intent(in) :: itype
    real*8, optional :: rho0(:), vel0(:)

    integer, parameter :: maxrad = 200

    real*8, allocatable :: rho(:), vel(:)
    integer :: inuc, jnuc, i, j, kk, k, il, ir
    integer :: nr, nang, nwav, l, l2, iy, jy, nangj, nrj
    real*8 :: h, rholm(maxrad), xx(maxrad), sum
    real*8 :: wang(302), rmid, v(3), rmax
    real*8 :: ylm(302,0:14,0:14)
    real*8 :: charge, bnd(10,-2:maxrad+3), r, q, dq, y, rfac
    real*8, allocatable :: rads(:,:), rq(:), rqq(:), wrads(:), d2r(:,:), trid2r(:,:)
    real*8, allocatable :: x(:,:,:)
    real*8 :: rs(302), xunit(3,302), dqs(302)
    integer :: ipvt(maxrad), ierr, i1, intqs(302), intq
    real*8, dimension(0:maxrad,0:14,0:14) :: ulm, a, b, c

    allocate(rho(mesh%n),vel(mesh%n))
    if (present(rho0)) then
       rho = rho0
    else
       if (itype == 0) then
          rho = mesh%rho(:,0)
       else
          call error('poiss_grid','unknown itype',2)
       end if
    end if
    mesh%vel = 0d0
    vel = 0d0

    ! fill rads
    allocate(rads(maxrad,mol%n),x(3,302,mol%n))
    allocate(wrads(maxrad),rq(maxrad),rqq(maxrad),d2r(maxrad,7),trid2r(maxrad,3))
    do inuc = 1, mol%n
       nr = z2nr(mol%z(inuc))
       nang = z2nang(mol%z(inuc))
       rmid = 1d0/real(mol%z(inuc),8)**third
       call rmesh(nr,rmid,rads(:,inuc),rq,rqq,wrads,d2r,trid2r)
       call ld0302(x(1,:,inuc),x(2,:,inuc),x(3,:,inuc),wang,nang)
    end do

    do inuc = 1, mol%n
       ! generate radial mesh
       nr = z2nr(mol%z(inuc))
       rmid = 1d0/real(mol%z(inuc),8)**third
       h = 1d0 / (nr+1)
       call rmesh(nr,rmid,rads(:,inuc),rq,rqq,wrads,d2r,trid2r)

       nang = z2nang(mol%z(inuc))
       nwav = z2lmax(mol%z(inuc))
       kk = mesh%ktab(inuc)
       charge = 0d0
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             charge = charge + mesh%w(kk)*rho(kk)
          end do
       end do

       call ld0302(x(1,:,inuc),x(2,:,inuc),x(3,:,inuc),wang,nang)
       call ycalc(x(:,:,inuc),nang,nwav,ylm(1:nang,0:nwav,0:nwav))

       do i = 0, nwav
          do j = 0, nwav
             kk = mesh%ktab(inuc)

             ! calculate the density rholm
             do ir = 1, nr
                sum = 0d0
                do il = 1, nang
                   sum = sum + wang(il) * ylm(il,i,j) * mesh%wnuc(kk+il) * rho(kk+il)
                end do
                rholm(ir) = sum
                kk = kk + nang
             end do
             
             ! solve the radial lm equation
             l = max(i,j)
             l2 = l*(l+1)
             
             ! lhs (linpack band matrix)
             do k = 1, 7
                do ir = 1, nr
                   bnd(11-k,ir+k-4)=d2r(ir,k)
                end do
             end do
             do ir = 1, nr
                xx(ir)=-fourpi*rads(ir,inuc)*rholm(ir)
                bnd(7,ir)=bnd(7,ir)-l2/rads(ir,inuc)**2
             end do
             
             ! infinite-r boundary condition
             if (l == 0) then
                xx(nr)=xx(nr)-charge*d2r(nr,5)
                xx(nr-1)=xx(nr-1)-charge*d2r(nr-1,6)
                xx(nr-2)=xx(nr-2)-charge*d2r(nr-2,7)
             end if
             call dgbfa(bnd(1,1),10,nr,3,3,ipvt,ierr)
             if(ierr.ne.0) call error('poiss_grid','error in dgbfa',2)
             call dgbsl(bnd(1,1),10,nr,3,3,ipvt,xx,0)

             ulm(0,i,j) = 0.d0
             do ir = 1, nr
                ulm(ir,i,j)=xx(ir)
             end do

             kk = mesh%ktab(inuc)
             do ir = 1, nr
                r = rads(ir,inuc)
                do il = 1, nang
                   vel(kk+il)=vel(kk+il)+ulm(ir,i,j)*ylm(il,i,j)/r
                end do
                kk = kk + nang
             end do

             if(l.eq.0)then
                call spline2(h,ulm(0,i,j),a(0,i,j),b(0,i,j),c(0,i,j),nr,2,0.d0,charge)
             else
                call spline2(h,ulm(0,i,j),a(0,i,j),b(0,i,j),c(0,i,j),nr,1,0.d0,0.d0)
             end if
          end do
       end do

       do jnuc = 1, mol%n
          if (inuc == jnuc) cycle
          nangj = z2nang(mol%z(jnuc))
          nrj = z2nr(mol%z(jnuc))

          ! interpolate only l=0 for this pair if the centers are too far apart
          v = mol%x(:,inuc) - mol%x(:,jnuc)
          r = sqrt(dot_product(v,v))
          rmax = 1.5d0 * (covrad(mol%z(inuc)) + covrad(mol%z(jnuc)))
          if (r > rmax) then
             nwav = 14
          else
             nwav = min(z2lmax(mol%z(jnuc)),14)
          end if

          ! need a lebedev run for all the atoms first!
          kk = mesh%ktab(jnuc)
          do ir = 1, nrj
             do il = 1, nangj
                v = mol%x(:,jnuc) + x(:,il,jnuc) * rads(ir,jnuc) - mol%x(:,inuc)
                rs(il) = sqrt(dot_product(v,v))
                xunit(:,il) = v / rs(il)
                q = rs(il) / (rs(il)+rmid)
                intqs(il) = int((nr+1)*q)
                dqs(il) = q-intqs(il)*h
             end do
             call ycalc(xunit,nangj,nwav,ylm(1:nangj,0:nwav,0:nwav)) ! note the nl discrepancy
             do iy=0,nwav
                do jy=0,nwav
                   do il=1,nangj
                      intq=intqs(il)
                      dq=dqs(il)
                      y=ulm(intq,iy,jy)
                      rfac=y+dq*(a(intq,iy,jy)+dq*(b(intq,iy,jy)+dq*c(intq,iy,jy)))
                      vel(kk+il)=vel(kk+il)+rfac*ylm(il,iy,jy)/rs(il)
                   end do
                end do
             end do
             kk=kk+nangj
          end do
       end do
    end do

    deallocate(wrads,rq,rqq,d2r,trid2r)
    deallocate(rho,rads,x)
    if (present(vel0)) then
       vel0 = vel
       deallocate(vel)
    else
       if (itype == 0) then
          if (allocated(mesh%vel)) deallocate(mesh%vel)
          call move_alloc(vel,mesh%vel)
       else
          call error('poiss_grid','unknown itype',2)
       end if
    endif

   end subroutine poiss_grid

end module meshmod
