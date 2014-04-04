! -*- mode: F90 -*-

!> Calculate the XDM three-body coefficients. This routine is
!> only valid for water trimers.
subroutine sandbox_xdmc9(mol,mesh)
  implicit none

  type(molecule), intent(in) :: mol !< Input molecule
  type(tmesh), intent(in) :: mesh !< Associated molecular mesh

  integer :: narg
  real*8 :: c1br, c2br
  character*(mline) :: line, hfword
  integer :: i, lp, ii, jj, kk
  real*8 :: ntotal
  integer :: j, k, k1, k2, istat, inuc
  real*8 :: d, atpol(mol%n), fac, rvdw, c6, c8, c10, rc
  real*8 :: c6com, c8com, c10com, xij(3), ifac, r, r1, r2
  real*8 :: e, f(3,mol%n), qfreq(3,mol%n,3,mol%n), qfac
  real*8 :: mm(3,mol%n), v(mol%n), q(mol%n)
  real*8 :: hirsh(mesh%n)
  real*8 :: xcm(3,3)
  real*8 :: rab, rac, rbc, ca, cb, cc, c9, ea, eb, ec

  narg = command_argument_count()
  call getarg(3,line)
  read (line,*,err=999,end=999) hfword
  lp = 1
  if (.not.isreal(chf,hfword,lp)) then
     if (trim(lower(hfword)) == "blyp") then
        chf = chf_blyp
     elseif (trim(lower(hfword)) == "b3lyp") then
        chf = chf_b3lyp
     elseif (trim(lower(hfword)) == "bhandhlyp" .or. trim(lower(hfword)) == "bhandh"&
        .or. trim(lower(hfword)) == "bhah" .or. trim(lower(hfword)) == "bhahlyp") then
        chf = chf_bhahlyp
     elseif (trim(lower(hfword)) == "camb3lyp" .or. trim(lower(hfword)) == "cam-b3lyp" ) then
        chf = chf_camb3lyp
     elseif (trim(lower(hfword)) == "pbe") then
        chf = chf_pbe
     elseif (trim(lower(hfword)) == "pbe0") then
        chf = chf_pbe0
     elseif (trim(lower(hfword)) == "lcwpbe" .or. trim(lower(hfword)) == "lc-wpbe") then
        chf = chf_lcwpbe
     elseif (trim(lower(hfword)) == "pw86" .or. trim(lower(hfword)) == "pw86pbe") then
        chf = chf_pw86
     elseif (trim(lower(hfword)) == "b971" .or. trim(lower(hfword)) == "b97-1") then
        chf = chf_b971
     else
        call error("sandbox_xdm","unknown functional",2)
     endif
  endif
  ntotal = sum(mesh%w * (mesh%rho(:,1)+mesh%rho(:,2)))
  if (abs(ntotal - mol%nelec) > 0.1d0) then
     write (iout,'("WARNING -- inconsistent nelec. I hope you know what you are doing.")')
  endif

  ! calculate the XDM-related properties
  mm = 0d0
  v = 0d0
  q = 0d0

  rewind(ihrsh)
  do inuc = 1, mol%n
     if (mol%z(inuc) < 1) cycle
     read(ihrsh) (hirsh(i),i=1,mesh%n)
     mm(:,inuc) = 0d0
     v(inuc) = 0d0
     q(inuc) = 0d0

     ! calculate hole dipole and moments
     do i = 1, mesh%n
        r = sqrt(dot_product(mesh%x(:,i)-mol%x(:,inuc),mesh%x(:,i)-mol%x(:,inuc)))
        r1 = max(0.d0,r-mesh%bxdm(i,1))
        r2 = max(0.d0,r-mesh%bxdm(i,2))

        mm(1,inuc)=mm(1,inuc) + mesh%w(i) * hirsh(i) * &
           (mesh%rho(i,1) * (r-r1)**2 + mesh%rho(i,2) * (r-r2)**2)
        mm(2,inuc)=mm(2,inuc) + mesh%w(i) * hirsh(i) * &
           (mesh%rho(i,1) * (r**2-r1**2)**2 + mesh%rho(i,2) * (r**2-r2**2)**2)
        mm(3,inuc)=mm(3,inuc) + mesh%w(i) * hirsh(i) * &
           (mesh%rho(i,1) * (r**3-r1**3)**2 + mesh%rho(i,2) * (r**3-r2**3)**2)
        v(inuc) = v(inuc) + mesh%w(i) * hirsh(i) * &
           (mesh%rho(i,1)+mesh%rho(i,2)) * r**3
        q(inuc) = q(inuc) + mesh%w(i) * hirsh(i) * (mesh%rho(i,1)+mesh%rho(i,2))
     enddo
  enddo

  ! calculate the dispersion energy from the moments
  do i = 1, mol%n
     if (mol%z(i) < 1) cycle
     atpol(i) = v(i) * frepol(mol%z(i)) / frevol(mol%z(i))
  enddo

  ! ONLY for water trimers !xxxx!

  ! list of c9 and geometric terms
  c9 = 0d0
  do i = 1, 3
     ea = 2d0/3d0 * mm(1,i) / atpol(i)
     do j = 4, 6
        eb = 2d0/3d0 * mm(1,j) / atpol(j)
        do k = 7, 9
           ec = 2d0/3d0 * mm(1,k) / atpol(k)
           c9 = 4d0/9d0 * mm(1,i)*mm(1,j)*mm(1,k) * (ea+eb+ec) / (ea+eb) / (ea+ec) / (eb+ec)
           rab = sqrt(dot_product(mol%x(:,j) - mol%x(:,i),mol%x(:,j) - mol%x(:,i)))
           rac = sqrt(dot_product(mol%x(:,k) - mol%x(:,i),mol%x(:,k) - mol%x(:,i)))
           rbc = sqrt(dot_product(mol%x(:,k) - mol%x(:,j),mol%x(:,k) - mol%x(:,j)))
           ca = dot_product(mol%x(:,j) - mol%x(:,i),mol%x(:,k) - mol%x(:,i)) / rab / rac
           cb = dot_product(mol%x(:,k) - mol%x(:,j),mol%x(:,i) - mol%x(:,j)) / rab / rbc
           cc = dot_product(mol%x(:,j) - mol%x(:,k),mol%x(:,i) - mol%x(:,k)) / rac / rbc
           fac = (3 * ca * cb * cc + 1) / rab**3 / rac**3 / rbc**3
           write (iout,'(1p,7(E15.8,X))') c9, fac, rab, rac, rbc
        end do
     end do
  end do

  !    ! center of mass for the three water moleculesa
  !    xcm = 0d0
  !    do i = 1, 3
  !       do j = 1, 3
  !          xcm(:,i) = xcm(:,i) + mol%x(:,3*(i-1)+j) * mol%z(3*(i-1)+j)
  !       end do
  !       xcm(:,i) = xcm(:,i) / 10d0
  !    end do
  !    rab = sqrt(dot_product(xcm(:,2) - xcm(:,1),xcm(:,2) - xcm(:,1)))
  !    rac = sqrt(dot_product(xcm(:,3) - xcm(:,1),xcm(:,3) - xcm(:,1)))
  !    rbc = sqrt(dot_product(xcm(:,3) - xcm(:,2),xcm(:,3) - xcm(:,2)))
  !    ca = dot_product(xcm(:,2) - xcm(:,1),xcm(:,3) - xcm(:,1)) / rab / rac
  !    cb = dot_product(xcm(:,3) - xcm(:,2),xcm(:,1) - xcm(:,2)) / rab / rbc
  !    cc = dot_product(xcm(:,2) - xcm(:,3),xcm(:,1) - xcm(:,3)) / rac / rbc
  !
  !    ! only for water trimers
  !    c9 = 0d0
  !    do ii = 1, 3
  !       i = ii
  !       ea = 2d0/3d0 * mm(1,i) / atpol(i)
  !       do jj = 1, 3
  !          j = 3 + jj
  !          eb = 2d0/3d0 * mm(1,j) / atpol(j)
  !          do kk = 1, 3
  !             k = 6 + kk
  !             ec = 2d0/3d0 * mm(1,k) / atpol(k)
  !
  !             c9 = c9 + 4d0/9d0 * mm(1,i)*mm(1,j)*mm(1,k) * &
  !                (ea+eb+ec) / (ea+eb) / (ea+ec) / (eb+ec)
  !
  !          end do
  !       end do
  !    end do
  !
  !    d = 0.52917720859d0
  !    fac = (3 * ca * cb * cc + 1) / rab**3 / rac**3 / rbc**3
  !    write (iout,'(1p,7(E15.8,X))') c9, fac, rab, rac, rbc

  return

999 write (iout,'(/"Error in command line: postg2 file.wfx xdmc9 functional"/)')
  stop 1

end subroutine sandbox_xdmc9
