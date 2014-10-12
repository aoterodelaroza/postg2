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

  character*(mline) :: line, msg, sint
  integer :: i, j, k
  real*8 :: ebr1, ux, b, alf, a, x, expx, enuc, eee, een
  real*8 :: etau, b1, alf1, a1, x1, ux1, v(3), r, exx, e, ec
  real*8 :: charge, zk, etotal, flip, r1, r2, zz, uu, fss, dsigs
  real*8 :: uxps1, uxps2, quads, m2, dn, m1
  real*8, allocatable :: vnuc(:), phi(:)
  logical :: mask(mprops)
  integer :: narg, iarg

  real*8, parameter :: tiny = 1d-20

  etotal = 0d0
  narg = command_argument_count()
  iarg = 2
  do while(.true.)
     iarg = iarg + 1
     if (iarg > narg) exit
     call getarg(iarg,line)
     line = adjustl(lower(line))
     if (trim(line) == "charge") then
        charge = sum(mesh%w * mesh%rho(:,0))
        write (iout,'("Charge (charge) = ",1p,E22.14)') charge
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
     elseif (trim(line) == "enn") then
        enuc = 0d0
        do i = 1, mol%n
           do j = i+1, mol%n
              v = mol%x(:,i) - mol%x(:,j)
              r = sqrt(dot_product(v,v))
              enuc = enuc + mol%z(i) * mol%z(j) / r
           end do
        end do
        write (iout,'("Nuclear repulsion energy (enn) = ",1p,E22.14)') enuc
        etotal = etotal + enuc
     elseif (trim(line) == "ekin") then
        etau = 0.5d0 * (sum(mesh%w * mesh%tau(:,1))+sum(mesh%w * mesh%tau(:,2)))
        write (iout,'("Kinetic energy (ekin) = ",1p,E22.14)') etau
        etotal = etotal + etau
     elseif (trim(line) == "ebr") then
        ebr1 = 0d0
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
                 ebr1 = ebr1 + mesh%w(i) * 0.5d0 * mesh%rho(i,j) * ux
              endif
           end do
        end do
        write (iout,'("Classic BR exchange (ebr) = ",1p,E22.14)') ebr1
        etotal = etotal + ebr1
     elseif (trim(line) == "eee") then
        mask = .false.
        mask(20) = .true.
        call propts_grid(mol,mesh,mask)
        eee = 0.5d0*sum(mesh%w * mesh%rho(:,0) * mesh%vel)
        write (iout,'("Coulomb energy (eee) = ",1p,E22.14)') eee
        etotal = etotal + eee
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
        een = sum(mesh%w * mesh%rho(:,0) * vnuc)
        deallocate(vnuc)
        write (iout,'("Coulomb energy (een) = ",1p,E22.14)') een
        etotal = etotal + een
     elseif (trim(line) == "exx") then
        mask = .false.
        mask(21) = .true.
        call propts_grid(mol,mesh,mask)
        exx = sum(mesh%w * mesh%exdens(:,0))
        write (iout,'("Exchange energy (exx) = ",1p,E22.14)') exx
        etotal = etotal + exx
     elseif (trim(line) == "ecb13_dyn_opp") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)
        ec = 0d0

        do k = 1, mesh%n
           flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
           flip = min(flip,1d0)
           uxps1 = -2d0*mesh%exdens(k,1)/max(mesh%rho(k,1),tiny)
           uxps2 = -2d0*mesh%exdens(k,2)/max(mesh%rho(k,2),tiny)
           r1 = mesh%xlns(k,1) / max(abs(uxps1),tiny)
           r2 = mesh%xlns(k,2) / max(abs(uxps1),tiny)
           zz = 0.63d0 * (r1 + r2)
           uu = -0.8d0 * mesh%rho(k,1) * mesh%rho(k,2) * zz**3 / (1d0+zz)
           ec = ec + mesh%w(k) * uu * (1d0-flip)
        end do

        write (iout,'("B13 dynamical correlation energy, opp (ecb13_dyn_opp) = ",1p,E22.14)') ec
        etotal = etotal + ec

     elseif (trim(line) == "ecb13_dyn_par") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        ec = 0d0
        do i = 1, 2
           do k = 1, mesh%n
              flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
              flip = min(flip,1d0)
              uxps1 = -2d0*mesh%exdens(k,i)/max(mesh%rho(k,i),tiny)
              r1 = mesh%xlns(k,i) / max(abs(uxps1),tiny)

              ! calculate the moments
              dsigs = mesh%tau(k,i) - 0.25d0 * mesh%drho2(k,i) / max(mesh%rho(k,i),tiny)
              ! quads = (mesh%d2rho(k,i)-2d0*dsigs) / 6d0
              ! call bhole(mesh%rho(k,i),quads,mesh%xlns(k,i),b,alf,a)
              ! x = alf * b
              m2 =  0.25d0 * mesh%xlns(k,i) / pi

              dn = 1 - mesh%xlns(k,i) - flip * mesh%xlns(k,mod(i,2)+1)
              fss = min(3d0*mesh%rho(k,i)*dn/(dsigs*m2),1d0)
              zz = 2d0 * 0.88d0 * r1
              ec = ec + mesh%w(k) * (1d0-fss) * mesh%rho(k,i) * dsigs * zz**5  / (1d0+0.5d0*zz)
           end do
        end do
        ec = -0.005d0 * ec
        write (iout,'("B13 dynamical correlation energy, parallel (ecb13_dyn_par) = ",1p,E22.14)') ec
        etotal = etotal + ec

     elseif (trim(line) == "ecb13_stat_opp") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        ec = 0d0
        do k = 1, mesh%n
           if (mesh%rho(k,1)<tiny .and. mesh%rho(k,2) < tiny) cycle
           flip = min((1d0-mesh%xlns(k,1))/max(mesh%xlns(k,2),tiny),(1d0-mesh%xlns(k,2))/max(mesh%xlns(k,1),tiny))
           flip = min(flip,1d0)
           uxps1 = -2d0*mesh%exdens(k,1)/max(mesh%rho(k,1),tiny)
           uxps2 = -2d0*mesh%exdens(k,2)/max(mesh%rho(k,2),tiny)
           ec = ec - 0.5d0 * mesh%w(k) * flip * (mesh%rho(k,1)*uxps2 + mesh%rho(k,2)*uxps1)
        end do
        write (iout,'("B13 static correlation energy, opp (ecb13_stat_opp) = ",1p,E22.14)') ec
        etotal = etotal + ec

     elseif (trim(line) == "ecb13_stat_par") then
        mask = .false.
        mask(1:16) = .true.
        mask(21) = .true.
        mask(22) = .true.
        call propts_grid(mol,mesh,mask)

        ec = 0d0
        do i = 1, 2
           do k = 1, mesh%n
              if (abs(mesh%rho(k,i)) < tiny) cycle
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
              fss = min(3d0*mesh%rho(k,i)*dn/(dsigs*m2),1d0)
              uu = -fss * dsigs / 3d0 / max(mesh%rho(k,i),tiny) * m1
              ec = ec + mesh%w(k) * mesh%rho(k,i) * uu
           end do
        end do
        ec = 0.5d0 * ec
        write (iout,'("B13 static correlation energy, parallel (ecb13_stat_par) = ",1p,E22.14)') ec
        etotal = etotal + ec

     elseif (trim(line) == "libxc") then
#ifdef HAVE_LIBXC
        iarg = iarg + 1
        call getarg(iarg,line)
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
