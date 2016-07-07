! -*- mode: F90 -*-
! Copyright (c) 2013-2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
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

!> Calculate the XDM dispersion energy and derivatives.
subroutine sandbox_xdm(mol,mesh)
  implicit none

  type(molecule), intent(in) :: mol !< Molecule
  type(tmesh), intent(in) :: mesh !< Associated molecular mesh.

  integer :: narg
  real*8 :: c1br, c2br
  character*(mline) :: line, hfword
  integer :: i, lp
  real*8 :: ntotal
  integer :: j, k1, k2, istat, inuc
  real*8 :: d, atpol(mol%n), fac, rvdw, c6, c8, c10, rc
  real*8 :: c6com, c8com, c10com, xij(3), ifac, r, r1, r2
  real*8 :: e, f(3,mol%n), qfreq(3,mol%n,3,mol%n), qfac
  real*8 :: mm(3,mol%n), v(mol%n), q(mol%n)
  real*8 :: hirsh(mesh%n), vtot
  real*8 :: sumc6

  call getarg(3,line)
  read (line,*,end=999,err=999) c1br
  call getarg(4,line)
  read (line,*,end=999,err=999) c2br

  write(IOUT,'("a1          ",F12.6)') c1br
  write(IOUT,'("a2(ang)     ",F12.6)') c2br
  c2br = c2br / 0.52917720859d0

  narg = command_argument_count()
  if (narg == 5) then
     call getarg(5,line)
     read (line,*,end=999,err=999) hfword
     lp = 1
     if (.not.isreal(chf,hfword,lp)) then
        if (trim(lower(hfword)) == "blyp") then
           chf = chf_blyp
           write(iout,'("a_hf        blyp")')
        elseif (trim(lower(hfword)) == "b3lyp") then
           chf = chf_b3lyp
           write(iout,'("a_hf        b3lyp")')
        elseif (trim(lower(hfword)) == "bhandhlyp" .or. trim(lower(hfword)) == "bhandh"&
           .or. trim(lower(hfword)) == "bhah" .or. trim(lower(hfword)) == "bhahlyp") then
           chf = chf_bhahlyp
           write(iout,'("a_hf        bhandhlyp")')
        elseif (trim(lower(hfword)) == "camb3lyp" .or. trim(lower(hfword)) == "cam-b3lyp" ) then
           chf = chf_camb3lyp
           write(iout,'("a_hf        camb3lyp")')
        elseif (trim(lower(hfword)) == "pbe") then
           chf = chf_pbe
           write(iout,'("a_hf        pbe")')
        elseif (trim(lower(hfword)) == "pbe0") then
           chf = chf_pbe0
           write(iout,'("a_hf        pbe0")')
        elseif (trim(lower(hfword)) == "lcwpbe" .or. trim(lower(hfword)) == "lc-wpbe") then
           chf = chf_lcwpbe
           write(iout,'("a_hf        lcwpbe")')
        elseif (trim(lower(hfword)) == "pw86" .or. trim(lower(hfword)) == "pw86pbe") then
           chf = chf_pw86
           write(iout,'("a_hf        pw86pbe")')
        elseif (trim(lower(hfword)) == "b971" .or. trim(lower(hfword)) == "b97-1") then
           chf = chf_b971
           write(iout,'("a_hf        b971")')
        else
           call error("sandbox_xdm","unknown functional",2)
        endif
     else
        write(iout,'("a_hf        ",f12.6)') chf
     endif
  else
     write(iout,'("a_hf        ",f12.6)') chf
  endif
  write (iout,'("mesh size     ",I10)') mesh%n
  write (iout,'("nelec         ",F12.6)') mol%nelec
  write (iout,'("nelec (promol)",F12.6)') sum(mesh%w * mesh%promol)
  write (iout,'("nelec, alpha  ",F12.6)') sum(mesh%w * mesh%rho(:,1))
  write (iout,'("nelec, beta   ",F12.6)') sum(mesh%w * mesh%rho(:,2))
  ntotal = sum(mesh%w * (mesh%rho(:,1)+mesh%rho(:,2)))
  write (iout,'("nelec, total  ",F12.6)') ntotal
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

  write (iout,'("moments and volumes ")')
  write (iout,'("# i At        <M1^2>             <M2^2>              <M3^2>           Volume              Vfree")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,5(E18.10,X))') i, ptable(mol%z(i)), mm(1,i), &
        mm(2,i), mm(3,i), v(i), frevol(mol%z(i))
  enddo
  write (iout,'("#")')

  write (iout,'("hirshfeld charges ")')
  write (iout,'("# i At        Charge")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,1(E18.10,X))') i, ptable(mol%z(i)), mol%z(i)-q(i)
  enddo
  write (iout,'("#")')

  ! calculate the dispersion energy from the moments
  do i = 1, mol%n
     if (mol%z(i) < 1) cycle
     atpol(i) = v(i) * frepol(mol%z(i)) / frevol(mol%z(i))
  enddo

  ! print out the polarizabilities
  write (iout,'("atomic polarizabilities ")')
  write (iout,'("# i At        Polariz.")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,F10.4)') i, ptable(mol%z(i)), atpol(i)
  enddo
  write (iout,'("#")')
  write (iout,'("molecular polarizability ",F12.6)') sum(atpol(1:mol%n))
  c6 = sum(atpol(:))*sum(mm(1,:))/2d0
  write (iout,'("molecular c6 ",F12.6)') c6
  write (iout,'("#")')

  ! calculate coefficients
  write (iout,'("coefficients and distances (a.u.)")')
  write (iout,'("# i  j       dij            C6               C8               C10              Rc           Rvdw")') 
  e = 0d0
  f = 0d0
  q = 0d0
  sumc6 = 0d0
  do i = 1, mol%n
     if (mol%z(i) < 1) cycle
     do j = i, mol%n
        if (mol%z(j) < 1) cycle
        xij = mol%x(:,j)-mol%x(:,i)
        d = sqrt(dot_product(xij,xij))
        fac = atpol(i)*atpol(j)/(mm(1,i)*atpol(j)+mm(1,j)*atpol(i))
        c6 = fac*mm(1,i)*mm(1,j)
        c8 = 1.5d0*fac*(mm(1,i)*mm(2,j)+mm(2,i)*mm(1,j))
        c10 = 2.d0*fac*(mm(1,i)*mm(3,j)+mm(3,i)*mm(1,j))&
           +4.2d0*fac*mm(2,i)*mm(2,j)
        rc = (sqrt(c8/c6) + sqrt(sqrt(c10/c6)) +&
           sqrt(c10/c8)) / 3.D0
        rvdw = c1br * rc + c2br
        if (i == j) then
           sumc6 = sumc6 + c6
        else
           sumc6 = sumc6 + 2d0 * c6
        endif
        if (d > 1d-5) then
           e = e - c6 / (rvdw**6 + d**6) - c8 / (rvdw**8+d**8) - &
              c10 / (rvdw**10 + d**10)
           c6com = 6.d0*c6*d**4/(rvdw**6+d**6)**2
           c8com = 8.d0*c8*d**6/(rvdw**8+d**8)**2
           c10com = 10.d0*c10*d**8/(rvdw**10+d**10)**2
           f(:,i) = f(:,i) + (c6com+c8com+c10com) * xij
           f(:,j) = f(:,j) - (c6com+c8com+c10com) * xij
           do k1 = 1, 3
              do k2 = 1, 3
                 if (k1 == k2) then
                    ifac = 1d0
                 else
                    ifac = 0d0
                 endif
                 qfac = &
                    c6com  * (-ifac - 4*xij(k1)*xij(k2)/d**2 + 12*xij(k1)*xij(k2)*d**4/(rvdw**6+d**6)) + &
                    c8com  * (-ifac - 6*xij(k1)*xij(k2)/d**2 + 16*xij(k1)*xij(k2)*d**6/(rvdw**8+d**8)) + &
                    c10com * (-ifac - 8*xij(k1)*xij(k2)/d**2 + 20*xij(k1)*xij(k2)*d**8/(rvdw**10+d**10)) 
                 qfreq(k1,i,k2,j) = qfac
                 qfreq(k2,j,k1,i) = qfac
              enddo
           enddo
        endif
        write (iout,'(I3,X,I3,1p,E14.6,X,3(E16.9,X),2(E13.6,X))') &
           i, j, d, c6, c8, c10, rc, rvdw
     end do
  end do
  write (iout,'("#")')
  write (iout,'("dimer c6 ",F12.6)') sumc6
  write (iout,'("#")')

  ! sum rules for the second derivatives
  do i = 1, mol%n
     do k1 = 1, 3
        do k2 = 1, 3
           qfreq(k1,i,k2,i) = 0d0
           do j = 1, mol%n
              if (j == i) cycle
              qfreq(k1,i,k2,i) = qfreq(k1,i,k2,i) - qfreq(k1,i,k2,j)
           enddo
        enddo
     enddo
  enddo

  write (iout,'("dispersion energy ",1p,E20.12)') e
  write (iout,'("scf energy ",1p,E20.12)') mol%escf
  write (iout,'("total energy (SCF+XDM) ",1p,E20.12)') mol%escf+e
  write (iout,'("dispersion forces ")')
  write (iout,'("# i          Fx                   Fy                   Fz")')
  do i = 1, mol%n
     write (iout,'(I3,X,1p,3(E20.12,X))') i, f(:,i)
  enddo
  write (iout,'("#")')
  write (iout,'("dispersion force constant matrix ")')
  write (iout,'("# i  xyz   j   xyz    Exixj ")')
  do i = 1, mol%n
     do k1 = 1, 3
        do j = 1, i-1
           do k2 = 1, 3
              write (iout,'(4(I3,X),1p,E20.12)') i, k1, j, k2, qfreq(k1,i,k2,j)
           enddo
        enddo
        do k2 = 1, k1
           write (iout,'(4(I3,X),1p,E20.12)') i, k1, i, k2, qfreq(k1,i,k2,j)
        enddo
     enddo
  enddo
  write (iout,'("#"/)')

  return

999 write (iout,'(/"Error in command line: postg2 file.wfx xdm a1 a2 functional"/)')
  stop 1

end subroutine sandbox_xdm
