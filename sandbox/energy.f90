! -*- mode: F90 -*-

!> Write the Hartree and the kinetic energy densities
!> to the standard output.
subroutine sandbox_energy(mol,mesh)
#ifdef HAVE_LIBXC
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
#endif
  use meshmod
  implicit none

  type(molecule), intent(in) :: mol !< Input geometry
  type(tmesh), intent(inout) :: mesh !< Molecular mesh

#ifdef HAVE_LIBXC
  type libxc_functional
     integer :: family ! LDA, GGA, etc.
     integer :: id     ! identifier
     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional
  type(libxc_functional) :: ifun
#endif

  character*40 :: carg(100)
  character*(mline) :: line, msg, sint
  integer :: i, j, k
  real*8 :: ux, b, alf, a, x, expx
  real*8 :: b1, alf1, a1, x1, ux1, v(3), r, e, rdum
  real*8 :: charge, zk, etotal, flip, r1, r2, zz, uu, fss, dsigs, zij
  real*8 :: uxps1, uxps2, quads, m2, dn, m1, ucstat, ucboth
  real*8 :: brmom1, brmom2, m3, m4, ucdii(2), ucdij, fii(2)
  real*8, allocatable :: vnuc(:), phi(:)
  logical :: mask(mprops), ok
  integer :: narg, iarg, idx, lp

  real*8, parameter :: tiny = 1d-20
  real*8, parameter :: small = 1d-10

  ! ! xxxx
  ! mask = .false.
  ! mask(1:22) = .true.
  ! call propts_grid(mol,mesh,mask)
  ! open(unit=35,name="bleh.dat")
  ! allocate(mesh%dsigs(1:mesh%n,0:2))
  ! allocate(mesh%quads(1:mesh%n,0:2))
  ! do i = 1, mesh%n
  !    read (35,'(999(E20.14,X))') mesh%w(i), mesh%rho(i,1), mesh%rho(i,2),&
  !       mesh%tau(i,1), mesh%tau(i,2), mesh%dsigs(i,1), mesh%dsigs(i,2),&
  !       mesh%quads(i,1), mesh%quads(i,2), mesh%exdens(i,1), mesh%exdens(i,2),&
  !       mesh%xlns(i,1), mesh%xlns(i,2)
  !    mesh%rho(i,0) = mesh%rho(i,1) + mesh%rho(i,2)
  !    mesh%tau(i,0) = mesh%tau(i,1) + mesh%tau(i,2)
  !    mesh%dsigs(i,0) = mesh%dsigs(i,1) + mesh%dsigs(i,2)
  !    mesh%quads(i,0) = mesh%quads(i,1) + mesh%quads(i,2)
  !    mesh%exdens(i,0) = mesh%exdens(i,1) + mesh%exdens(i,2)
  ! end do
  ! close(35)

  ! build the argument array
  carg = "  "
  narg = 0
  do iarg = 3, command_argument_count()
     call getarg(iarg,line)
     narg = narg + 1
     carg(narg) = trim(adjustl(line))
  end do

  etotal = 0d0
  iarg = 0
  do while(.true.)
     iarg = iarg + 1
     if (iarg > narg) exit
     line = lower(carg(iarg))
     if (trim(line) == "charge") then
        charge = sum(mesh%w * mesh%rho(:,0))
        write (iout,'("Charge (charge) = ",1p,E22.14)') charge
        e = 0d0
     elseif (trim(line) == "mocharge") then
        allocate(phi(mesh%n))
        rewind(iimo)
        do i = 1, mol%nmo
           read(iimo) phi(1:mesh%n)
           charge = sum(mesh%w * phi**2)
           write (msg,'(I99)') i
           write (iout,'("MO (mocharge, ",A,") = ",1p,E22.14)') trim(adjustl(msg)), charge
        end do
        deallocate(phi)
        e = 0d0
     elseif (trim(line) == "enn") then
        e = 0d0
        do i = 1, mol%n
           do j = i+1, mol%n
              v = mol%x(:,i) - mol%x(:,j)
              r = sqrt(dot_product(v,v))
              e = e + mol%z(i) * mol%z(j) / r
           end do
        end do
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Nuclear repulsion energy (enn) = ",1p,E22.14)') e 
        etotal = etotal + e 
     elseif (trim(line) == "ekin") then
        e = 0.5d0 * (sum(mesh%w * mesh%tau(:,1))+sum(mesh%w * mesh%tau(:,2)))
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Kinetic energy (ekin) = ",1p,E22.14)') e
        etotal = etotal + e 
     elseif (trim(line) == "ebr") then
        e = 0d0
        do j = 1,2
           do i = 1, mesh%n
              b = mesh%bxdm(i,j)
              b1 = b
              alf = mesh%alf(i,j)
              alf1 = alf
              a = mesh%prefac(i,j)
              a1 = a
              x = alf * b
              x1 = x
              expx = exp(-x)
              if (b > 1d-12) then
                 ux = -1d0/b * (1 - expx - 0.5d0*x*expx)
                 ux1 = ux
                 e = e + mesh%w(i) * 0.5d0 * mesh%rho(i,j) * ux
              endif
           end do
        end do
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Classic BR exchange (ebr) = ",1p,E22.14)') e 
        etotal = etotal + e
     elseif (trim(line) == "eee") then
        mask = .false.
        mask(20) = .true.
        call propts_grid(mol,mesh,mask)
        e = 0.5d0*sum(mesh%w * mesh%rho(:,0) * mesh%vel)
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Coulomb energy (eee) = ",1p,E22.14)') e
        etotal = etotal + e
     elseif (trim(line) == "een") then
        allocate(vnuc(mesh%n))
        vnuc = 0d0
        do i = 1, mol%n
           do j = 1, mesh%n
              v = mesh%x(:,j) - mol%x(:,i)
              r = sqrt(dot_product(v,v))
              vnuc(j) = vnuc(j) - mol%z(i)/r
           end do
        end do
        e = sum(mesh%w * mesh%rho(:,0) * vnuc)
        deallocate(vnuc)
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Coulomb energy (een) = ",1p,E22.14)') e
        etotal = etotal + e
     elseif (trim(line) == "exx") then
        mask = .false.
        mask(21) = .true.
        call propts_grid(mol,mesh,mask)
        e = sum(mesh%w * mesh%exdens(:,0))
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("Exchange energy (exx) = ",1p,E22.14)') e
        etotal = etotal + e
     elseif (trim(line) == "ecb13_dyn_opp") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)
        e = 0d0

        do k = 1, mesh%n
           if (mesh%rho(k,1)<small .and. mesh%rho(k,2) < small) cycle
           flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
           flip = min(flip,1d0)
           uxps1 = -2d0*mesh%exdens(k,1)/max(mesh%rho(k,1),tiny)
           uxps2 = -2d0*mesh%exdens(k,2)/max(mesh%rho(k,2),tiny)
           r1 = mesh%xlns(k,1) / max(abs(uxps1),tiny)
           r2 = mesh%xlns(k,2) / max(abs(uxps1),tiny)
           zz = 0.63d0 * (r1 + r2)
           uu = -0.8d0 * mesh%rho(k,1) * mesh%rho(k,2) * zz**3 / (1d0+zz)
           e = e + mesh%w(k) * uu * (1d0-flip)
        end do

        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("B13 dynamical correlation energy, opp (ecb13_dyn_opp) = ",1p,E22.14)') e
        etotal = etotal + e

     elseif (trim(line) == "ecb13_dyn_par") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        e = 0d0
        do i = 1, 2
           do k = 1, mesh%n
              if (mesh%rho(k,i)<small) cycle
              flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
              flip = min(flip,1d0)
              uxps1 = -2d0*mesh%exdens(k,i)/max(mesh%rho(k,i),tiny)
              r1 = mesh%xlns(k,i) / max(abs(uxps1),tiny)

              ! calculate the moments
              dsigs = mesh%tau(k,i) - 0.25d0 * mesh%drho2(k,i) / max(mesh%rho(k,i),tiny)
              m2 =  0.25d0 * mesh%xlns(k,i) / pi

              dn = 1 - mesh%xlns(k,i) - flip * mesh%xlns(k,mod(i,2)+1)
              fss = min(3d0*mesh%rho(k,i)*dn/max(dsigs*m2,tiny),1d0)
              zz = 2d0 * 0.88d0 * r1
              e = e + mesh%w(k) * (1d0-fss) * mesh%rho(k,i) * dsigs * zz**5  / (1d0+0.5d0*zz)
           end do
        end do
        e = -0.005d0 * e
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("B13 dynamical correlation energy, parallel (ecb13_dyn_par) = ",1p,E22.14)') e
        etotal = etotal + e

     elseif (trim(line) == "ecb13_stat_opp") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        e = 0d0
        do k = 1, mesh%n
           if (mesh%rho(k,1)<small .and. mesh%rho(k,2) < small) cycle
           flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
           flip = min(flip,1d0)
           uxps1 = -2d0*mesh%exdens(k,1)/max(mesh%rho(k,1),tiny)
           uxps2 = -2d0*mesh%exdens(k,2)/max(mesh%rho(k,2),tiny)
           e = e - 0.5d0 * mesh%w(k) * flip * (mesh%rho(k,1)*uxps2 + mesh%rho(k,2)*uxps1)
        end do
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("B13 static correlation energy, opp (ecb13_stat_opp) = ",1p,E22.14)') e
        etotal = etotal + e

     elseif (trim(line) == "ecb13_stat_par") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        e = 0d0
        do i = 1, 2
           do k = 1, mesh%n
              if (mesh%rho(k,i)<small) cycle
              flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
              flip = min(flip,1d0)
              uxps1 = -2d0*mesh%exdens(k,i)/max(mesh%rho(k,i),tiny)
              r1 = mesh%xlns(k,i) / max(abs(uxps1),tiny)

              ! calculate the moments
              dsigs = mesh%tau(k,i) - 0.25d0 * mesh%drho2(k,i) / max(mesh%rho(k,i),tiny)
              quads = (mesh%d2rho(k,i)-2d0*dsigs) / 6d0
              call bhole(mesh%rho(k,i),quads,mesh%xlns(k,i),b,alf,a)
              x = alf * b
              m1 = mesh%xlns(k,i) * (0.25d0-exp(-x)*(0.125d0*x+0.25d0)) / b / pi
              m2 = 0.25d0 * mesh%xlns(k,i) / pi

              dn = 1 - mesh%xlns(k,i) - flip * mesh%xlns(k,mod(i,2)+1)
              fss = min(3d0*mesh%rho(k,i)*dn/max(dsigs*m2,tiny),1d0)
              uu = -fss * dsigs / 3d0 / max(mesh%rho(k,i),tiny) * m1
              e = e + mesh%w(k) * mesh%rho(k,i) * uu
           end do
        end do
        e = 0.5d0 * e
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (iout,'("B13 static correlation energy, parallel (ecb13_stat_par) = ",1p,E22.14)') e
        etotal = etotal + e

     elseif (trim(line) == "ecb13_xtra2" .or. trim(line) == "ecb13_xtra3") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)
        if (trim(line) == "ecb13_xtra2") then 
           idx = 2
        else
           idx = 3
        endif

        e = 0d0
        do k = 1, mesh%n
           ucstat = 0d0
           ucboth = 0d0
           if (mesh%rho(k,1)<small .and. mesh%rho(k,2) < small) cycle

           ! parallel spins
           flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
           flip = min(flip,1d0)
           ucdii = 0d0
           zij = 0d0
           do i = 1, 2
              dsigs = mesh%tau(k,i) - 0.25d0 * mesh%drho2(k,i) / max(mesh%rho(k,i),tiny)
              quads = (mesh%d2rho(k,i)-2d0*dsigs) / 6d0
              uxps1 = -2d0*mesh%exdens(k,i)/max(mesh%rho(k,i),tiny)
              r1 = mesh%xlns(k,i) / max(abs(uxps1),tiny)
              zij = zij + r1
              m2 =  0.25d0 * mesh%xlns(k,i) / pi
              dn = 1 - mesh%xlns(k,i) - flip * mesh%xlns(k,mod(i,2)+1)
              call bhole(mesh%rho(k,i),quads,mesh%xlns(k,i),b,alf,a)
              x = alf * b
              m3 = 0.25d0 + 1d0 / x**2 
              m3 = m3 - exp(-x) * (0.25d0+1d0/x) / x
              m3 = mesh%xlns(k,i) * b * m3 / pi
              brmom1 = 4d0*pi*m3
              m4 = 0.25d0 + 3.0d0 / x**2
              m4 = m4 * b**2 / pi
              m4 = mesh%xlns(k,i) * m4
              brmom2 = 4d0 * pi * m4
              fii(i) = min(3.d0*mesh%rho(k,i)*dn/(dsigs*brmom2),1.d0)
              ucstat = ucstat - fii(i) * dsigs * brmom1/6.D0
              zz = 2d0 * 0.88d0 * r1
              ucdii(i) = -0.005d0 * mesh%rho(k,i) * dsigs * zz**5/(1.d0+0.5d0*zz)
           end do
           zij = 0.63d0 * zij 

           ! opposite spins
           uxps1 = -2d0*mesh%exdens(k,1)/max(mesh%rho(k,1),tiny)
           uxps2 = -2d0*mesh%exdens(k,2)/max(mesh%rho(k,2),tiny)
           ucstat = ucstat - 0.5d0*mesh%rho(k,1)*flip*uxps2 - 0.5d0*mesh%rho(k,2)*flip*uxps1
           ucdij = -0.8d0 * mesh%rho(k,1) * mesh%rho(k,2) * zij**3/(1.d0+zij)
           ucboth = ucstat + (1d0-flip) * ucdij + (1d0-fii(1)) * ucdii(1) + (1d0-fii(2)) * ucdii(2)

           if (ucstat > -small) cycle
           x = ucstat / ucboth
           e = e + mesh%w(k) * x**idx * ucboth
        end do
        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        if (trim(line) == "ecb13_xtra2") then 
           write (iout,'("B13 xtra2 = ",1p,E22.14)') e
        else
           write (iout,'("B13 xtra3 = ",1p,E22.14)') e
        end if
        etotal = etotal + e

     elseif (trim(line) == "b13") then
        carg(iarg+13:iarg+24) = carg(iarg+1:iarg+12)
        carg(iarg+1 ) = "enn"
        carg(iarg+2 ) = "1.0"
        carg(iarg+3 ) = "ekin"
        carg(iarg+4 ) = "1.0"
        carg(iarg+5 ) = "eee"
        carg(iarg+6 ) = "1.0"
        carg(iarg+7 ) = "een"
        carg(iarg+8 ) = "1.0"
        carg(iarg+9 ) = "exx"
        carg(iarg+10) = "1.0"
        carg(iarg+11) = "ecb13_stat_opp"
        carg(iarg+12) = "0.55242281838985"
        carg(iarg+13) = "ecb13_stat_par"
        carg(iarg+14) = "0.84448814464331"
        carg(iarg+15) = "ecb13_dyn_opp"
        carg(iarg+16) = "0.64038922588770"
        carg(iarg+17) = "ecb13_dyn_par"
        carg(iarg+18) = "0.55917730862225"
        carg(iarg+19) = "ecb13_xtra2"
        carg(iarg+20) = "0.82462761325139"
        carg(iarg+21) = "ecb13_xtra3"
        carg(iarg+22) = "-0.37974035284994"
        narg = iarg + 22
     elseif (trim(line) == "libxc") then
#ifdef HAVE_LIBXC
        iarg = iarg + 1
        line = carg(iarg)
        read (line,*) ifun%id
        ifun%family = xc_f90_family_from_id(ifun%id)
        if (ifun%family < 0) call error('sandbox_energy','invalid functional',2)
        if (ifun%id == 21) cycle ! exchange in 1D
        if (ifun%id == 160 .or. ifun%id == 182) cycle ! van Leeuwen & Baerends not implemented
        if (ifun%id == 207)  cycle ! BJ06 modified potential
        if (ifun%id == 208)  cycle ! Tran-Blaha BJ06
        if (ifun%id == 209)  cycle ! Rasanen, Pittali, Proetto BJ06

        call xc_f90_func_init(ifun%conf,ifun%info,ifun%id,XC_UNPOLARIZED)
        if (ifun%id == XC_LDA_C_XALPHA) &
           call xc_f90_lda_c_xalpha_set_par(ifun%conf,0.d0)
        call xc_f90_info_name(ifun%info,msg)

        if (mol%wfntyp == 1) call error('energy','spin-polarized libxc not implemented',2)

        e = 0d0
        do i = 1, mesh%n
           if (mesh%rho(i,0) < tiny) cycle
           select case(ifun%family)
           case (XC_FAMILY_LDA)
              call xc_f90_lda_exc(ifun%conf, 1, mesh%rho(i,0), zk)
           case (XC_FAMILY_GGA)
              call xc_f90_gga_exc(ifun%conf, 1, mesh%rho(i,0), mesh%drho2(i,0), zk)
           case (XC_FAMILY_MGGA)
              call xc_f90_mgga_exc(ifun%conf, 1, mesh%rho(i,0), mesh%drho2(i,0),&
                 mesh%d2rho(i,0), 0.5d0 * mesh%tau(i,0), zk)
           end select
           e = e + mesh%w(i) * mesh%rho(i,0) * zk
        end do

        call xc_f90_func_end(ifun%conf)

        lp = 1
        ok = isreal(rdum,carg(iarg+1),lp)
        if (.not.ok) then
           rdum = 1d0
        else
           iarg = iarg + 1
        endif
        e = e * rdum
        write (sint,'(I99)') ifun%id
        write (iout,'("Exc: ",A," (libxc ",A,") = ",1p,E22.14)') &
           trim(adjustl(msg)), trim(adjustl(sint)), e
#else
        call error("sandbox_energy","postg2 not compiled with libxc",2)
#endif
        etotal = etotal + e
     else
        call error("sandbox_energy","unknown keyword",2)
     endif
  end do
  if (etotal /= 0d0) then
     write (iout,'("Total energy (etotal) = ",1p,E22.14)') etotal
  endif

end subroutine sandbox_energy

