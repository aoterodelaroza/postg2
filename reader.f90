! Copyright (c) 2013-2015 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
! and Axel D. Becke <axel.becke@dal.ca>
!
! postg is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module reader
  implicit none

  public
  public :: reader_init ! initialize the matrices
  public :: readwfn ! create a molecule from a wfn file
  public :: readwfx ! create a molecule from a wfx file
  public :: readfchk ! create a molecule from a fchk file
  public :: readmolden ! create a molecule from a molden file
  public :: readtck ! create a molecule from a terachem checkpoint file

  real*8, private :: m5d6d(6,5)
  real*8, private :: m7d10d(10,7)

  ! gaussian and molden seq: XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
  integer, parameter :: d10gaumap(11:20) = (/11, 12, 13, 17, 14, 15, 18, 19, 16, 20/)

contains

  subroutine reader_init()

    m5d6d(:,1) = (/        -0.5d0,         -0.5d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    m5d6d(:,2) = (/         0.0d0,          0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
    m5d6d(:,3) = (/         0.0d0,          0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /)
    m5d6d(:,4) = (/ sqrt(3d0)/2d0, -sqrt(3d0)/2d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    m5d6d(:,5) = (/         0.0d0,          0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 /)

    m7d10d(:,1) = (/0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-1.5d0/sqrt(5d0),0.0d0,0.0d0,-1.5d0/sqrt(5d0),0.0d0/)
    m7d10d(:,2) = (/-sqrt(3d0/8d0),0.0d0,0.0d0,-sqrt(3d0/40d0),0.0d0,0.0d0,sqrt(6d0/5d0),0.0d0,0.0d0,0.0d0/)
    m7d10d(:,3) = (/0.0d0,-sqrt(3d0/8d0),0.0d0,0.0d0,-sqrt(3d0/40d0),0.0d0,0.0d0,sqrt(6d0/5d0),0.0d0,0.0d0/)
    m7d10d(:,4) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,sqrt(3d0)/2d0,0.0d0,0.0d0,-sqrt(3d0)/2d0,0.0d0/)
    m7d10d(:,5) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/)
    m7d10d(:,6) = (/sqrt(5d0/8d0),0.0d0,0.0d0,-3d0/sqrt(8d0),0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)
    m7d10d(:,7) = (/0.0d0,-sqrt(5d0/8d0),0.0d0,0.0d0,3d0/sqrt(8d0),0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)
  
  end subroutine reader_init

  !> Read a gaussian wfn file
  function readwfn(file) result(m)
    use types
    use param

    character*(mline), intent(in) :: file !< Name of the wfn file
    type(molecule) :: m !< Resulting molecule

    character*4 :: orbtyp
    character*2 :: dums, elem
    integer :: i, j, istat, imax, icount, ioc, num1, num2
    real*8 :: zreal, ene, ene0
    logical :: isfrac
    integer :: nalpha
    character*8 :: dum1, dum3, dum4
    character*18 :: dum2
    character*1024 :: line

    integer :: idum, idum2
    integer, parameter :: luwfn = 10

    ! set title
    m%name = file
    m%useecp = .false.

    ! read number of atoms, primitives, orbitals
    open(luwfn,file=file,status='old')
    read (luwfn,*)
    read (luwfn,101) orbtyp, m%nmo, m%npri, m%n

    ! atomic positions and numbers
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for atomic positions',2)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for atomic numbers',2)
    m%x = 0d0
    m%z = 0
    m%charge = 0d0
    do i = 1, m%n
       read(luwfn,106) elem, m%x(:,i), zreal
       m%charge = m%charge + zreal
       m%z(i) = elem2z(elem)
       if (m%z(i) /= nint(zreal)) m%useecp = .true.
    end do

    ! center assignments, types of primitives
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for itype',2)
    read(luwfn,102) (m%icenter(i),i=1,m%npri)
    read(luwfn,102) (m%itype(i),i=1,m%npri)
    if (any(m%itype(1:m%npri) > 56)) then
       call error("readwfn","primitive type not supported",2)
    endif

    ! primitive exponents
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for exponents',2)
    read(luwfn,103) (m%e(i),i=1,m%npri)

    ! deal with ecps
    dums=""
    do while (dums.ne."MO")
       read (luwfn,'(A2)') dums
    enddo
    backspace(luwfn)

    ! occupations and orbital coefficients
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for occupations',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory&
       & for orbital coefficients',2)
    m%mult = 1
    isfrac = .false.
    num1 = 0
    num2 = 0
    ene0 = -1d30
    nalpha = -1
    m%nelec = 0d0
    do i = 1, m%nmo
       read(luwfn,104) m%occ(i), ene
       read(luwfn,105) (m%c(i,j),j=1,m%npri)
       m%nelec = m%nelec + m%occ(i)
       ioc = nint(m%occ(i))
       if (abs(ioc-m%occ(i)) > 1d-10) then
          isfrac = .true.
       else if (ioc == 1) then
          num1 = num1 + 1
       else if (ioc == 2) then
          num2 = num2 + 1
       endif
       if (ene < ene0-1d-3) nalpha = i-1
       ene0 = ene
    end do
    read(luwfn,*) dum1
    read(luwfn,'(A80)') line
    line = trim(adjustl(line))
    line = line(index(line,'=')+1:)
    read (line,*) m%escf

    ! figure out charge and multiplicity
    ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
    m%charge = m%charge - m%nelec
    if (isfrac) then
       m%wfntyp = 3
       call error("readwfn","natural orbital wfn files not supported"&
          &,2)
    else if (num1 == 0) then
       m%wfntyp = 0
    else if (num2 == 0) then
       m%wfntyp = 1
       if (nalpha == -1) then
          m%mult = nint(m%nelec) + 1
       else
          m%mult = nalpha - (nint(m%nelec) - nalpha) + 1
       endif
       if (m%mult < 0) call error("readwfn","nbeta > nalpha",2)
    else
       m%wfntyp = 2
       m%mult = count(nint(m%occ(1:m%nmo)) == 1) + 1
    endif
    close(luwfn)

101 format (4X,A4,10X,3(I5,15X))
102 format(20X,20I3)
103 format(10X,5E14.7)
104 format(35X,F12.7,15X,F12.6)
105 format(5(E16.8))
106 format(2X,A2,20X,3F12.8,10X,F5.1)

  end function readwfn

  !> Read a gaussian wfx file
  function readwfx(file) result(m)
    use types
    use param
    use io

    character*(mline), intent(in) :: file !< Name of the wfx file
    type(molecule) :: m !< Resulting molecule

    integer :: i, j, istat, ncore, kk, lp, idum
    real*8 :: zreal
    character*(mline) :: line, tag
    logical :: keyw(8)

    integer, parameter :: luwfn = 10

    ! set title
    m%name = file
    m%useecp = .false.

    ! first pass
    open (luwfn,file=file,status='old')
    m%n = 0
    m%nmo = 0
    m%charge = 0
    m%mult = 0
    ncore = 0
    m%npri = 0
    m%escf = 0d0
    do while (.true.)
       read(luwfn,'(A)',end=10) line
       line = adjustl(line)
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Nuclei>") then
             read (luwfn,*) m%n
          elseif (trim(line) == "<Number of Occupied Molecular Orbitals>") then
             read (luwfn,*) m%nmo
          elseif (trim(line) == "<Net Charge>") then
             read (luwfn,*) m%charge
          elseif (trim(line) == "<Electronic Spin Multiplicity>") then
             read (luwfn,*) m%mult
          elseif (trim(line) == "<Number of Core Electrons>") then
             read (luwfn,*) ncore
          elseif (trim(line) == "<Number of Primitives>") then
             read(luwfn,*) m%npri
          elseif (trim(line) == "<Energy = T + Vne + Vee + Vnn>") then
             read(luwfn,*) m%escf
          endif
       endif
    enddo
10  continue

    if (m%n == 0) call error("readwfx","Number of Nuclei tag not found",2)
    if (m%nmo == 0) call error("readwfx","Number of Occupied Molecular Orbitals tag not found",2)
    if (m%mult == 0) call error("readwfx","Electronic Spin Multiplicity tag not found",2)
    if (m%npri == 0) call error("readwfx","Number of Primitives tag not found",2)
    if (ncore > 0) m%useecp = .true.

    ! allocate memory
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for atomic positions',2)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for atomic numbers',2)
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for exponents',2)
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for occupations',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for orbital coefficients',2)
    allocate(m%eps(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for orbital energies',2)

    ! second pass
    rewind(luwfn)
    keyw = .false.
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line = adjustl(line)
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Atomic Numbers>") then
             m%z = read_integers(luwfn,m%n)
             keyw(1) = .true.
          elseif (trim(line) == "<Nuclear Cartesian Coordinates>") then
             m%x = reshape(read_reals1(luwfn,3*m%n),shape(m%x))
             keyw(2) = .true.
          elseif (trim(line) == "<Primitive Centers>") then
             m%icenter = read_integers(luwfn,m%npri)
             keyw(3) = .true.
          elseif (trim(line) == "<Primitive Types>") then
             m%itype = read_integers(luwfn,m%npri)
             if (any(m%itype(1:m%npri) > 56)) &
                call error("readwfx","primitive type not supported",2)
             keyw(4) = .true.
          elseif (trim(line) == "<Primitive Exponents>") then
             m%e = read_reals1(luwfn,m%npri)
             keyw(5) = .true.
          elseif (trim(line) == "<Molecular Orbital Occupation Numbers>") then
             m%occ = read_reals1(luwfn,m%nmo)
             m%nelec = sum(m%occ)
             keyw(6) = .true.
          elseif (trim(line) == "<Molecular Orbital Energies>") then
             m%eps = read_reals1(luwfn,m%nmo)
             keyw(7) = .true.
          elseif (trim(line) == "<Molecular Orbital Primitive Coefficients>") then
             read(luwfn,*)
             do i = 1, m%nmo
                read(luwfn,*)
                read(luwfn,*)
                m%c(i,:) = read_reals1(luwfn,m%npri)
             enddo
             keyw(8) = .true.
          endif
       endif
    enddo
20  continue
    if (any(.not.keyw)) call error("readwfx","missing array in wfx file",2)

    ! wavefuntion type
    ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
    if (m%mult == 1) then
       m%wfntyp = 0
    else
       if (any(m%occ > 1)) then
          m%wfntyp = 2
       else
          m%wfntyp = 1
       endif
    end if

    close(luwfn)

  end function readwfx

  !> Read a gaussian fchk file
  function readfchk(file) result(m)
    use types
    use param

    character*(mline), intent(in) :: file !< Name of the wfn file
    type(molecule) :: m !< Resulting molecule

    character*(mline) :: line
    integer :: lp, idum, nalpha, nbeta, nshel, ncshel, lmax, nbas
    integer :: istat, i, j, k, l, ifac
    integer :: acent, nn, nm, nl, jj
    logical :: ok, isbeta, isecp
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:)
    real*8, allocatable :: xat(:), exppri(:), ccontr(:), pccontr(:), mocoef(:)
    real*8 :: acoef
    real*8, allocatable :: dcoef(:,:), scoef(:,:)

    integer, parameter :: luwfn = 10

    ! set title
    m%name = file
    m%useecp = .false.

    ! first pass: dimensions
    isbeta = .false.
    isecp = .false.
    open(luwfn,file=file,status='old')
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line = adjustl(line)
       lp = 45
       if (line(1:6) == "Charge") then
          ok = isinteger(idum,line,lp)
          m%charge = idum
       elseif (line(1:12) == "Multiplicity") then
          ok = isinteger(m%mult,line,lp)
       elseif (line(1:19) == "Number of electrons") then
          ok = isinteger(idum,line,lp)
          m%nelec = idum
       elseif (line(1:15) == "Number of atoms") then
          ok = isinteger(m%n,line,lp)
       elseif (line(1:25) == "Number of alpha electrons") then
          ok = isinteger(nalpha,line,lp)
       elseif (line(1:25) == "Number of basis functions") then
          ok = isinteger(nbas,line,lp)
       elseif (line(1:24) == "Number of beta electrons") then
          ok = isinteger(nbeta,line,lp)
       elseif (line(1:27) == "Number of contracted shells") then
          ok = isinteger(ncshel,line,lp)
       elseif (line(1:26) == "Number of primitive shells") then
          ok = isinteger(nshel,line,lp)
       elseif (line(1:24) == "Highest angular momentum") then
          ok = isinteger(lmax,line,lp)
       elseif (line(1:12) == "Total Energy") then
          ok = isreal(m%escf,line,lp)
       elseif (line(1:21) == "Beta Orbital Energies") then
          isbeta = .true.
       elseif (line(1:8) == "ECP-LMax") then
          isecp = .true.
       endif
    enddo
20  continue

    if (.not.isbeta) then
       m%wfntyp = 0
    else
       m%wfntyp = 1
    endif
    if (isecp) call error("readfchk","ECPs not supported.",2)

    ! Count the number of MOs
    if (m%wfntyp == 0) then
       m%nmo = m%nelec / 2
       allocate(m%occ(m%nmo),stat=istat)
       if (istat /= 0) call error('readfchk','could not allocate memory for occ',2)
       m%occ = 2d0
    else if (m%wfntyp == 1) then
       m%nmo = m%nelec
       allocate(m%occ(m%nmo),stat=istat)
       if (istat /= 0) call error('readfchk','could not allocate memory for occ',2)
       m%occ = 1d0
    endif

    ! second pass
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for shell types',2)
    allocate(m%x(3,m%n),m%z(m%n),xat(3*m%n),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for geometry',2)
    allocate(exppri(nshel),ccontr(nshel),pccontr(nshel),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for prim. shells',2)
    allocate(mocoef(nbas*m%nmo),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for MO coefs',2)
    rewind(luwfn)
    do while (.true.)
       read(luwfn,'(A)',end=30) line
       line = adjustl(line)
       lp = 45
       if (line(1:11) == "Shell types") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlt(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Number of primitives per shell") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlpri(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Shell to atom map") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlat(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:29) == "Current cartesian coordinates") then
          do i = 0, (3*m%n-1)/5
             read(luwfn,'(5E16.8)',end=30) (xat(5*i+j),j=1,min(5,3*m%n-5*i))
          enddo
       elseif (line(1:14) == "Atomic numbers") then
          do i = 0, (m%n-1)/6
             read(luwfn,'(6I12)',end=30) (m%z(6*i+j),j=1,min(6,m%n-6*i))
          enddo
       elseif (line(1:19) == "Primitive exponents") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (exppri(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:24) == "Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (ccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:31) == "P(S=P) Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (pccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:21) == "Alpha MO coefficients") then
          do i = 0, (nalpha*nbas-1)/5
             read(luwfn,'(5E16.8)',end=30) (mocoef(5*i+j),j=1,min(5,nalpha*nbas-5*i))
          enddo
       elseif (line(1:21) == "Beta MO coefficients") then
          do i = 0, (nbeta*nbas-1)/5
             read(luwfn,'(5E16.8)',end=30) (mocoef(nalpha*nbas+5*i+j),j=1,min(5,nbeta*nbas-5*i))
          enddo
       endif
    enddo
30  continue

    if (any(abs(ishlt) > 3)) &
       call error("readfchk","primitives > f not supported",2)

    ! geometry
    m%x = reshape(xat,shape(m%x))

    ! Count the number of primitives
    m%npri = 0
    do i = 1, ncshel
       if (ishlt(i) == 0) then
          ifac = 1
       else if (ishlt(i) == 1) then
          ifac = 3
       else if (ishlt(i) == -1) then
          ifac = 4
       else if (abs(ishlt(i)) == 2) then
          ifac = 6
       else if (abs(ishlt(i)) == 3) then
          ifac = 10
       endif
       m%npri = m%npri + ifac * ishlpri(i)
    enddo

    ! Assign primitive center and type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for coeffs',2)
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       acent = ishlat(i)
       if (ishlt(i) == 0) then
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             m%icenter(nn) = acent
             m%itype(nn) = 1
             m%e(nn) = exppri(nm+k)
             do l = 1, m%nmo
                m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * mocoef((l-1)*nbas+nl)
             end do
          end do
       else if (ishlt(i) == 1) then
          do j = 2, 4
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == -1) then
          do j = 1, 4
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                if (j /= 1) then
                   acoef = pccontr(nm+k)
                else
                   acoef = ccontr(nm+k)
                endif
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * acoef * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == 2) then
          do j = 5, 10
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == -2) then
          allocate(dcoef(5,m%nmo),scoef(6,m%nmo))
          do l = 1, m%nmo
             do j = 1, 5
                dcoef(j,l) = mocoef((l-1)*nbas+nl+j)
             end do
          end do
          scoef = matmul(m5d6d,dcoef)

          do j = 5, 10
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * scoef(j-4,l)
                end do
             end do
          end do
          deallocate(dcoef)
          deallocate(scoef)
          nl = nl + 5
       else if (ishlt(i) == 3) then
          do j = 11, 20
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == -3) then
          allocate(dcoef(7,m%nmo),scoef(10,m%nmo))
          do l = 1, m%nmo
             do j = 1, 7
                dcoef(j,l) = mocoef((l-1)*nbas+nl+j)
             end do
          end do
          scoef = matmul(m7d10d,dcoef)

          do j = 11, 20
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = d10gaumap(j)
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k) * scoef(j-10,l)
                end do
             end do
          end do
          deallocate(dcoef,scoef)
          nl = nl + 7
       endif
       nm = nm + ishlpri(i)
    end do

    deallocate(ishlt,ishlpri,ishlat)
    deallocate(xat)
    deallocate(exppri,ccontr,pccontr)
    deallocate(mocoef)
    close(luwfn)

  end function readfchk

  !> Read molden file
  function readmolden(file) result(m)
    use types
    use param

    character*(mline), intent(in) :: file
    type(molecule) :: m

    integer, parameter :: luwfn = 10

    character*(mline) :: line, word, word1, word2, keyword, wrest
    integer :: nbas, ncshel, nshel
    integer :: i, j, k, l, idum, idum1, ni, nj
    real*8 :: rdum, norm
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:)
    real*8, allocatable :: exppri(:), ccontr(:), mocoef(:)
    integer :: istat, nm, ifac, itypmax, isend
    logical :: isang
    integer :: acent, nn, nl, nalpha
    logical :: ok, isalpha
    logical :: is5d, is7f, is9g
    real*8, allocatable :: dcoef(:,:), scoef(:,:)

    line = ""
    m%name = file
    m%useecp = .false.
    is5d = .false.
    is7f = .false.
    is9g = .false.

    ! parse the molden file, first pass
    open(luwfn,file=file,status='old')

    do while(next_keyword())
       if (trim(lower(keyword)) == "atoms") then
          ! read the geometry
          m%n = 0
          read(luwfn,'(A)',end=40) line
          do while(index(lower(line),"[") == 0 .and. len(trim(line)) > 0)
             m%n = m%n + 1
             read(luwfn,'(A)',end=40) line
          end do
       elseif (trim(lower(keyword)) == "5d") then
          is5d = .true.
          is7f = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(lower(keyword)) == "5d7f") then
          is5d = .true.
          is7f = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(lower(keyword)) == "5d10f") then
          is5d = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(lower(keyword)) == "7f") then
          is7f = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(lower(keyword)) == "9g") then
          is9g = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(lower(keyword)) == "gto") then
          ! read the basis set details
          ncshel = 0
          nshel = 0
          do i = 1, m%n
             read(luwfn,'(A)',end=40) line
             read(luwfn,'(A)',end=40) line
             do while (index(line,".") /= 0)
                read (line,*) word, idum, rdum
                ncshel = ncshel + 1
                nshel = nshel + idum
                do j = 1, idum
                   read(luwfn,'(A)',end=40) line
                end do
                read(luwfn,'(A)',end=40) line
             end do
          end do
       elseif (trim(lower(keyword)) == "mo") then
          ! read the number of molecular orbitals
          nbas = 0
          m%nmo = 0
          m%nelec = 0
          nalpha = 0
          do while(.true.)
             read(luwfn,'(A)',end=30) line
             if (index(lower(line),"[") /= 0) exit
             if (index(lower(line),"ene=") /= 0) then
                ! spin
                read(luwfn,'(A)',end=30) line
                read (line,*) word1, word2
                isalpha = (trim(lower(word2)) == "alpha")
                ! occup
                read(luwfn,'(A)',end=30) line
                read (line,*) word, idum
                if (idum == 1) then
                   m%wfntyp = 1
                   m%nmo = m%nmo + 1
                   m%nelec = m%nelec + idum
                   if (isalpha) then
                      nalpha = nalpha + 1
                   endif
                elseif (idum == 2) then
                   m%wfntyp = 0
                   m%nmo = m%nmo + 1
                   m%nelec = m%nelec + idum
                elseif (idum == 0) then
                   continue
                else
                   call error('readmolden','wrong occupation',2)
                endif
                if (nbas == 0) then
                   read(luwfn,'(A)',end=30) line
                   do while (index(line,".") /= 0)
                      read (line,*) idum
                      read(luwfn,'(A)',end=30) line
                   end do
                   nbas = idum
                end if
             end if
          end do
30        continue
       else
          read(luwfn,'(A)',end=40) line
       end if
    end do

    ! allocate stuff
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for occ',2)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for shell types',2)
    allocate(exppri(nshel),ccontr(nshel),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for prim. shells',2)
    allocate(mocoef(nbas*m%nmo),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for MO coefs',2)

    ! fix the occupations for now -- i'll clean up later
    if (m%wfntyp == 0) then
       m%occ = 2
    elseif (m%wfntyp == 1) then
       m%occ = 1
    endif

    ! rewind
    rewind(luwfn)

    ! read the geometry
    read(luwfn,'(A)',end=40) line
    do while(.not.trim(lower(line)) == "[molden format]")
       read(luwfn,'(A)',end=40) line
    end do

    ! geometry header
    read(luwfn,'(A)',end=40) line
    read(line,*) word1, word2
    do while(.not.trim(lower(word1)) == "[atoms]")
       read(luwfn,'(A)',end=40) line
       read(line,*) word1, word2
    end do
    isang = (trim(lower(word2)) == "(angs)".or.trim(lower(word2)) == "(ang)")

    ! the actual geometry
    read(luwfn,'(A)',end=40) line
    allocate(m%x(3,m%n),m%z(m%n),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for geometry',2)
    do i = 1, m%n
       read(line,*) word1, idum1, m%z(i), m%x(:,i)
       read(luwfn,'(A)',end=40) line
    end do
    if (isang) m%x = m%x / 0.52917720859d0

    ! calculate charge and multiplicity... no energy info
    m%charge = 0
    do i = 1, m%n
       m%charge = m%charge + m%z(i)
    end do
    m%charge = m%charge - m%nelec
    m%mult = abs(nalpha - (m%nelec - nalpha)) + 1

    ! basis set
    do while(.not.trim(lower(line)) == "[gto]")
       read(luwfn,'(A)',end=40) line
    end do
    ni = 0
    nj = 0
    do i = 1, m%n
       read(luwfn,'(A)',end=40) line
       read(luwfn,'(A)',end=40) line
       do while (index(line,".") /= 0)
          ni = ni + 1
          read(line,*) word, idum
          do j = 1, idum
             nj = nj + 1
             read(luwfn,'(A)',end=40) line
             read(line,*) exppri(nj), ccontr(nj)
          end do
          ishlat(ni) = i
          ishlpri(ni) = idum
          if (lower(trim(word)) == "s") then
             ishlt(ni) = 0
          else if (lower(trim(word)) == "p") then
             ishlt(ni) = 1
          else if (lower(trim(word)) == "sp") then
             call error("readmolden","can't handle SP in gamess format",2)
          else if (lower(trim(word)) == "d") then
             if (is5d) then
                ishlt(ni) = -2
             else
                ishlt(ni) = 2
             endif
          else if (lower(trim(word)) == "f") then
             if (is7f) then
                ishlt(ni) = -3
             else
                ishlt(ni) = 3
             endif
          else
             call error("readmolden","basis set type not supported",2)
          endif
          read(luwfn,'(A)',end=40) line
       end do
    end do

    ! Count the number of primitives
    m%npri = 0
    do i = 1, ncshel
       if (ishlt(i) == 0) then
          ifac = 1
       else if (ishlt(i) == 1) then
          ifac = 3
       else if (ishlt(i) == -1) then
          ifac = 4
          call error('readmolden','ishlt = -1 not supported',2)
       else if (abs(ishlt(i)) == 2) then
          ifac = 6
       else if (abs(ishlt(i)) == 3) then
          ifac = 10
       endif
       m%npri = m%npri + ifac * ishlpri(i)
    enddo

    ! advance to the MO coefficients
    do while(index(lower(line),"[mo]") == 0)
       read(luwfn,'(A)',end=40) line
    end do
    do i = 1, m%nmo
       do while(.true.)
          read(luwfn,'(A)',end=40) line
          if (index(lower(line),"occup=") /= 0) then
             read(line,*) word, idum
             if (idum > 0) then
                do j = 1, nbas
                   read(luwfn,*,end=40) idum, mocoef((i-1)*nbas+j)
                end do
                exit
             end if
          end if
       end do
    end do

    ! Assign primitive center and type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for coeffs',2)

    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       acent = ishlat(i)
       if (ishlt(i) == 0) then
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             m%icenter(nn) = acent
             m%itype(nn) = 1
             m%e(nn) = exppri(nm+k)
             do l = 1, m%nmo
                m%c(l,nn) =  ccontr(nm+k) * mocoef((l-1)*nbas+nl)
             end do
          end do
       else if (ishlt(i) == 1) then
          do j = 2, 4
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == 2) then
          do j = 5, 10
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == -2) then
          allocate(dcoef(5,m%nmo),scoef(6,m%nmo))
          do l = 1, m%nmo
             do j = 1, 5
                dcoef(j,l) = mocoef((l-1)*nbas+nl+j)
             end do
          end do
          scoef = matmul(m5d6d,dcoef)

          do j = 5, 10
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = ccontr(nm+k) * scoef(j-4,l)
                end do
             end do
          end do
          deallocate(dcoef,scoef)
          nl = nl + 5 
       else if (ishlt(i) == 3) then
          do j = 11, 20
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = d10gaumap(j)
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = ccontr(nm+k) * mocoef((l-1)*nbas+nl)
                end do
             end do
          end do
       else if (ishlt(i) == -3) then
          allocate(dcoef(7,m%nmo),scoef(10,m%nmo))
          do l = 1, m%nmo
             do j = 1, 7
                dcoef(j,l) = mocoef((l-1)*nbas+nl+j)
             end do
          end do
          scoef = matmul(m7d10d,dcoef)

          do j = 11, 20
             do k = 1, ishlpri(i)
                nn = nn + 1
                m%icenter(nn) = acent
                m%itype(nn) = d10gaumap(j)
                m%e(nn) = exppri(nm+k)
                do l = 1, m%nmo
                   m%c(l,nn) = ccontr(nm+k) * scoef(j-10,l)
                end do
             end do
          end do
          deallocate(dcoef,scoef)
          nl = nl + 7
       endif
       nm = nm + ishlpri(i)
    end do

!    do j = 1, nn
!       write (*,'(3(I3,X),1p,2E20.12)') j, m%icenter(j), &
!          m%itype(j), m%e(j), m%c(1,j)
!    end do
!    stop 1

    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef)

    close(luwfn)

    return
40  continue
    call error("readmolden","unexpected end of file",2)

  contains

    function next_keyword()

      integer :: istart, iend, io
      logical :: next_keyword

      keyword = ""
      wrest = ""
      next_keyword = .false.

      do while(index(lower(line),"[") == 0)
         read(luwfn,'(A)',iostat=io) line
         if (io /= 0) return
      end do
      next_keyword =.true.
      istart = index(lower(line),"[") + 1
      iend = index(lower(line),"]") - 1
      keyword = line(istart:iend)
      wrest = line(iend+2:)

    end function next_keyword

  end function readmolden

  !> Read a terachem fchk file
  function readtck(file) result(m)
    use types
    use param

    character*(mline), intent(in) :: file !< Name of the wfn file
    type(molecule) :: m !< Resulting molecule

    character*(mline) :: line, word, word2
    integer :: lp, idum, nalpha, nbeta, nshel, ncshel, lmax, nbas
    integer :: istat, i, j, k, l, ifac, ii, jj, itypmax
    integer :: acent, nn, nm, nl
    logical :: ok
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:), npribas(:), nbaspri(:)
    real*8, allocatable :: xat(:), exppri(:), ccontr(:), pccontr(:), mocoef(:)
    real*8 :: aexp, acoef, rdum, norm

    integer, parameter :: luwfn = 10

    ! set title
    m%name = file
    m%useecp = .false.

    ! parse the tck file, first pass
    open(luwfn,file=file,status='old')
    read(luwfn,*,end=30) m%n
    read(luwfn,*,end=30) idum
    m%nelec = idum
    read(luwfn,*,end=30) idum
    m%charge = idum
    read(luwfn,*,end=30) m%mult
    read(luwfn,*,end=30) nbas
    read(luwfn,*,end=30) m%nmo
    read(luwfn,*,end=30) ncshel
    read(luwfn,*,end=30) nshel
    read(luwfn,'(/)',end=30)

    ! wfntyp -> closed-shell
    m%wfntyp = 0
    if (m%wfntyp == 0) then
       allocate(m%occ(m%nmo),stat=istat)
       if (istat /= 0) call error('readtck','could not allocate memory for occ',2)
       m%occ = 2d0
    endif

    ! allocates
    allocate(m%x(3,m%n),m%z(m%n),xat(3*m%n),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for geometry',2)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for shell types',2)
    allocate(exppri(nshel),ccontr(nshel),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for prim. shells',2)
    allocate(mocoef(nbas*m%nmo),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for MO coefs',2)

    ! geometry
    do i = 1, m%n
       read(luwfn,*,end=30) word, rdum, m%x(:,i)
       m%z(i) = nint(rdum)
    enddo

    read(luwfn,'(/)',end=30)
    do i = 1,m%n
       read(luwfn,*,end=30)
       nn = 0
       do while(.true.)
          read(luwfn,'(A80)',end=30) word
          lp = 1
          ok = isinteger(idum,word,lp)
          if (ok) then
             nn = nn + 1
             read(word,*,end=30) ii, word2, jj, aexp, acoef
             if (trim(word2) == "S") then
                ishlt(ii) = 0
             else if (trim(word2) == "P") then
                ishlt(ii) = 1
             else if (trim(word2) == "SP") then
                call error("readtck","can't handle SP in gamess format",2)
             else if (trim(word2) == "D") then
                ishlt(ii) = 2
             else if (trim(word2) == "F") then
                ishlt(ii) = 3
             else
                call error("readtck","basis set type not supported",2)
             endif
             ishlat(ii) = i
             exppri(jj) = aexp
             ccontr(jj) = acoef
          else
             if (nn > 0) ishlpri(ii) = nn
             nn = 0
             if (trim(word) /= "") then
                exit
             endif
          endif
       enddo
    enddo

    ! energy
    read (luwfn,*) m%escf

    ! advance to the MO
    word = ""
    do while(word /= "Occupied")
       read(luwfn,*,end=30) word
    end do
    do i = 1, nbas
       read(luwfn,*,end=30) (mocoef(j*nbas+i),j=0,m%nmo-1)
    enddo

    ! Count the number of primitives
    m%npri = 0
    do i = 1, ncshel
       if (ishlt(i) == 0) then
          ifac = 1
       else if (ishlt(i) == 1) then
          ifac = 3
       else if (ishlt(i) == -1) then
          ifac = 4
          call error('readtck','ishlt = -1 not supported',2)
       else if (ishlt(i) == 2) then
          ifac = 6
       else if (ishlt(i) == 3) then
          ifac = 10
       endif
       m%npri = m%npri + ifac * ishlpri(i)
    enddo

    ! Normalize the primitive shells
    nm = 0
    do i = 1, ncshel
       l = ishlt(i)
       norm = 0d0
       do j = 1, ishlpri(i)
          do k = 1, ishlpri(i)
             norm = norm + ccontr(nm+j) * ccontr(nm+k) * &
                sqrt((2*exppri(nm+j))**(l+1.5d0)) * sqrt((2*exppri(nm+k))**(l+1.5d0)) / &
                (exppri(nm+j) + exppri(nm+k))**(l+1.5d0)
          end do
       end do
       norm = sqrt(norm)
       do j = 1, ishlpri(i)
          ccontr(nm+j) = ccontr(nm+j) / norm
       end do
       nm = nm + ishlpri(i)
    end do

    ! Assign primitive center and type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for coeffs',2)
    allocate(npribas(m%npri),nbaspri(nbas),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for pri/bas index',2)
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       acent = ishlat(i)
       if (ishlt(i) == 0) then
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             npribas(nn) = nl
             nbaspri(nl) = nn
             m%icenter(nn) = acent
             m%itype(nn) = 1
             m%e(nn) = exppri(nm+k)
             m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
          end do
       else if (ishlt(i) == 1) then
          do j = 2, 4
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       else if (ishlt(i) == 2) then
          do j = 5, 10
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                ! terachem sequence is xy,xz,yz,xx,yy,zz
                if (j <= 7) then
                   m%itype(nn) = j+3
                else
                   m%itype(nn) = j-3
                endif
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       else if (ishlt(i) == 3) then
          do j = 11, 20
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       endif
       nm = nm + ishlpri(i)
    end do

    ! re-order MO coefficients in input
    itypmax = maxval(m%itype)
    nl = 0
    do l = 0, 3
       do i = 1, nbas
          ok = (l==0) .and. (m%itype(nbaspri(i))==1)
          ok = ok .or. (l==1) .and. (m%itype(nbaspri(i))>=2.and.m%itype(nbaspri(i))<=4)
          ok = ok .or. (l==2) .and. (m%itype(nbaspri(i))>=5.and.m%itype(nbaspri(i))<=10)
          ok = ok .or. (l==3) .and. (m%itype(nbaspri(i))>=11.and.m%itype(nbaspri(i))<=20)
          if (ok) then
             nl = nl + 1
             do k = 1, m%npri
                if (npribas(k) /= i) cycle
                do j = 1,m%nmo
                   m%c(j,k) = m%c(j,k) * mocoef((j-1)*nbas+nl)
                end do
             end do
          endif
          nn = nn + npribas(i)
       end do
    enddo

    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef,npribas,nbaspri)
    close(luwfn)

    return
30  continue
    call error("readtck","unexpected end of file",2)

  end function readtck

  !> Normalization constant for input primitives
  function gnorm(type,a) result(N)
    use tools_math
    use param
    integer, intent(in) :: type !< Type of primitive
    real*8, intent(in) :: a !< Primitive exponent
    real*8 :: N !< Normalization constant

    if (type == 1) then
       N = 2**(3d0/4d0) * a**(3d0/4d0) / pi**(3d0/4d0)
    else if (type >= 2 .and. type <= 4) then
       N = 2**(7d0/4d0) * a**(5d0/4d0) / pi**(3d0/4d0)
    else if (type >= 5 .and. type <= 7) then
       ! 5  6  7  8  9  10
       ! XX,YY,ZZ,XY,XZ,YZ
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type >= 7 .and. type <= 10) then
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0)
    else if (type >= 11 .and. type <= 13) then
       !  11  12  13  14  15  16  17  18  19  20
       ! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
       N = (2d0*a/pi)**0.75d0 * sqrt((8*a)**3 * fac(3) / fac(6))
!       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(15d0)
!       call error("gnorm","fixme: f primitives",2)
    else if (type >= 14 .and. type <= 19) then
       N = (2d0*a/pi)**0.75d0 * sqrt((8*a)**3 * fac(2) / fac(4) / fac(2))
!       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
!       call error("gnorm","fixme: f primitives",2)
    else if (type == 20) then
       N = (2d0*a/pi)**0.75d0 * sqrt((8*a)**3 * 1d0 / fac(2)**3)
!       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0)
!       call error("gnorm","fixme: f primitives",2)
    else
       call error("gnorm","fixme: primitive type not supported",2)
    endif

  endfunction gnorm

end module reader
