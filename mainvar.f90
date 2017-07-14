
       module mainvar


!------linkages.

!      used by - [program] main


!-----remarks.

!     Contains a list of variables and related functions for main


       implicit none


!------computational numbers.

       integer, parameter :: dl = kind(1.d0)


!------numerical constants.

       real(dl), parameter :: pi = 3.1415926535897932_dl
       real(dl), parameter :: boltzmann = 8.6173324e+04_dl !eV/(10^9 K)


!------physical constants.

       real(dl), parameter :: mpl = 1.220932e+22_dl !Planck mass in MeV
       real(dl), parameter :: amu = 931.494061_dl !atomic mass unit in MeV
       real(dl), parameter :: hbar = 6.582119514e-22_dl !Planck constant over 2\pi in MeV*s
       real(dl), parameter :: zeta3 = 1.2020569031_dl !Riemann zeta function of argument 3
       real(dl) :: deltamnp = 1.29333217_dl !mass difference in MeV between proton and neutron
       real(dl), parameter :: xmelec = 0.510998928_dl !electron mass in MeV
       real(dl) :: alphafs = 7.2973525698e-03_dl !fine structure constant


!------input values.

       real(dl) temp0k !Temperature at current epoch in 10^9 K


!------flags for output files.

       logical bbnflag !flag for bbn.dat
       logical equilflag !flag for equil.dat
       logical yeflag !flag for ye.dat
       integer runind
       logical saveflag

       logical mswflag !flag to use msw
       real(dl) maxepslinsave !maximum value of epsilon in eps_bins

       logical transflag !if .true., calculates neutrino transport

       logical filefound !used for finding existing files

       real(dl), dimension(:), allocatable :: xi !neutrino degeneracy parameters.
       integer nnu !number of total neutrinos


       type nuvar

         real(dl) :: dilfact
         real(dl) :: mass
         real(dl) :: mix !\sin(2\theta)
         real(dl), dimension(:), pointer :: occfrac

       end type nuvar

       type rhovar

         real(dl) :: rhoc0 !Criticial energy density at current epoch
         real(dl) :: og !Contribution from photons
         real(dl) :: om !matter
         real(dl) :: ol !dark energy
         real(dl) :: orad !extra radn. energy
         real(dl) :: temp0ev !plasma temp. at current epoch in eV
         real(dl) :: on !nu's, only use if onflag = .true.
         logical :: onflag
         type (nuvar), dimension(:,:), pointer :: nuparm

       end type rhovar

       type bin_scheme

         real(dl), dimension(:), pointer :: abscissas
         real(dl), dimension(:), pointer :: weights

       end type bin_scheme

       integer nbins !number of bins in nu spectrum
       type(bin_scheme) eps_bins !for neutrinos
       real(dl) och2 !oc*h**2, cold dark matter
       real(dl) h !hubble parameter at current epoch


!------Boole's Rule quantities.

       real(dl), dimension(4) :: brcoeffs = (/64._dl,24._dl,64._dl,28._dl/)/45._dl !coefficients
       real(dl), dimension(2) :: brends = (/14._dl,14._dl/)/45._dl !end-point coefficients


!------gamma function quantities.

       real(dl), dimension(6) :: gamcoeffs = (/76.18009172947146_dl &
                                             ,-86.50532032941677_dl &
                                             , 24.01409824083091_dl &
                                             ,-1.231739572450155_dl &
                                             , 1.208650973866179e-03_dl&
                                             ,-5.395239384953e-06_dl/)
       real(dl) :: stp = 2.5066282746310005_dl

       character(len=2) :: frepostr

!------MPI variables.

#ifdef prllel
       integer rank
       integer nprocs
#endif


       contains


!--------------------------------------------------------------------

       real(dl) function hubcalc(a,hparm)


!------linkages.

!      called by - [subroutine] integ_a
!                  [function] xe_de_h, xe_de_he
!      calls     - occfraclagint


!------remarks.

!      Calculates the Hubble rate using the Friedmann equation.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       type(rhovar), intent(in) :: hparm !recomvar variable


!------local variables.

       real(dl) rhoc0 !Criticial energy density at current epoch
       real(dl) og !Contribution from photons
       real(dl) om !matter
       real(dl) ol !dark energy
       real(dl) orad !extra radn. energy
       real(dl) on !ultra-relativistic \nu's, only use if onflag = .true.
       real(dl) tdil !diluted temperature factor
       real(dl) temp0ev !plasma temp in eV at current epoch
       real(dl) tempnuev !nu temp in eV at scale factor a
       real(dl) rhonu !total nu energy density in MeV
       real(dl) ms

       integer m,n


!------procedure.

       rhoc0 = hparm%rhoc0
       og = hparm%og
       om = hparm%om
       ol = hparm%ol
       orad  = hparm%orad

       if (hparm%onflag) then
         on = hparm%on
         hubcalc = sqrt(8._dl*pi/3._dl/mpl**2*rhoc0 &
                        *((og + orad + on)/a**4 + om/a**3 + ol))
       else
         temp0ev = hparm%temp0ev
         rhonu = 0._dl
         do n=1,2
           do m=1,size(hparm%nuparm,1)
             !tempnuev = temp0ev*(4._dl/11._dl)**(1._dl/3._dl)/a
             tempnuev = temp0ev*hparm%nuparm(m,n)%dilfact/a
             ms = hparm%nuparm(m,n)%mass/tempnuev
             rhonu = rhonu + 1._dl/2._dl/pi**2*tempnuev**4*1.e-24_dl &
                             *occfracint(hparm%nuparm(m,n)%occfrac, ms)
           end do
         end do
         hubcalc = sqrt(8._dl*pi/3._dl/mpl**2*(rhoc0 &
                        *((og + orad)/a**4 + om/a**3 + ol) + rhonu))
       end if

     
       return
       

       end function hubcalc


!--------------------------------------------------------------------

       real(dl) function dhubda(a,hparm)


!------linkages.

!      called by - [function] tmat_de
!      calls     - drholagint


!------remarks.

!      Calculates the derivative of the Hubble rate.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       type(rhovar), intent(in) :: hparm !recomvar variable


!------local variables.

       real(dl) rhoc0 !Criticial energy density at current epoch
       real(dl) og !Contribution from photons
       real(dl) om !matter
       real(dl) ol !dark energy
       real(dl) orad !extra radn. energy
       real(dl) on !ultra-relativistic nu's, only use if onflag = .true.
       real(dl) tdil !diluted temperature factor
       real(dl) temp0ev !plasma temp in eV at current epoch
       real(dl) tempnuev !nu temp in eV at scale factor a
       real(dl) rhonu !total nu energy density in MeV
       real(dl) ms
       real(dl) hub !part of the Hubble rate
       real(dl) drhoda !derivative of rhonu

       integer m,n


!------procedure.

       rhoc0 = hparm%rhoc0
       og = hparm%og
       om = hparm%om
       ol = hparm%ol
       orad  = hparm%orad

       if (hparm%onflag) then
         on = hparm%on
         hub = (og + orad + on)/a**4 + om/a**3 + ol
         dhubda = sqrt(8._dl*pi/3._dl/mpl**2*rhoc0)/2._dl/sqrt(hub) &
                  *(-4._dl*(og + orad + on)/a**5 - 3._dl*om/a**4)
       else
         temp0ev = hparm%temp0ev
         rhonu = 0._dl
         drhoda = 0._dl
         do n=1,2
           do m=1,size(hparm%nuparm,1)
             !tempnuev = temp0ev*(4._dl/11._dl)**(1._dl/3._dl)/a
             tempnuev = temp0ev*hparm%nuparm(m,n)%dilfact/a
             ms = hparm%nuparm(m,n)%mass/tempnuev
             rhonu = rhonu + 1._dl/2._dl/pi**2*tempnuev**4*1.e-24_dl &
                             *occfracint(hparm%nuparm(m,n)%occfrac, ms)
             drhoda = drhoda + 1._dl/2._dl/pi**2*tempnuev**4*1.e-24_dl &
                               *drhoint(hparm%nuparm(m,n)%occfrac, ms)/a
           end do
         end do
         drhoda = drhoda - 4._dl*rhonu/a
         hub = rhoc0*((og + orad)/a**4 + om/a**3 + ol) + rhonu
         dhubda = sqrt(8._dl*pi/3._dl/mpl**2)/2._dl/sqrt(hub)* &
                        (rhoc0*(-4._dl*(og + orad)/a**5 - 3._dl*om/a**4) + drhoda)
       end if

     
       return
       

       end function dhubda


!--------------------------------------------------------------------

       real(dl) function occfracint(moccfrac, ms)


!------linkages.

!      called by - [function] hubcalc
!      calls     - none


!------remarks.

!      Uses Gauss-Laguerre integration to compute nu energy density given
!      occupation fractions.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: moccfrac !occupation fractions in mass basis
       real(dl), intent(in) :: ms !dimensionless mass


!------local variables.

       integer k
       real(dl) integral, x, xint


!------procedure.

       integral = 0._dl
       do k=1,size(moccfrac) !integral over momentum
         x = eps_bins%abscissas(k)
         xint = x**2*sqrt(x**2 + ms**2)*moccfrac(k)
         integral = integral + xint*eps_bins%weights(k)
       end do

       occfracint = integral

     
       return
       

       end function occfracint


!--------------------------------------------------------------------

       real(dl) function drhoint(moccfrac, ms)


!------linkages.

!      called by - [function] dhubda
!      calls     - none


!------remarks.

!      Uses generic integration scheme to compute derivative of nu energy density given
!      occupation fractions.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: moccfrac !occupation fractions in mass basis
       real(dl), intent(in) :: ms !dimensionless mass


!------local variables.

       integer k
       real(dl) integral, x, xint


!------procedure.

       integral = 0._dl
       do k=1,size(moccfrac) !integral over momentum
         x = eps_bins%abscissas(k)
         if (x.ne.0._dl) then
           xint = 1._dl/sqrt(x**2 + ms**2)*x**2*moccfrac(k)*ms**2
         else
           xint = 0._dl
         end if
         integral = integral + xint*eps_bins%weights(k)
       end do

       drhoint = integral

     
       return
       

       end function drhoint


!-------------------------------------------------------------------     
       
       real(dl) function be_equil_calc(eps,tratio,etadeg)


!------linkages.

!      called by - [subroutine] bbn_no_therm_v2
!      calls     - none


!------remarks.

!      Calculates Bose-Einstein occupation probability.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: eps !input epsilon quantity
       real(dl), intent(in) :: tratio !input temperature ratio
       real(dl), intent(in) :: etadeg !input degeneracy parameter


!------procedure.


       if ((eps.gt.0._dl).and.(eps.le.300._dl)) then
         be_equil_calc = 1._dl/(exp(eps*tratio - etadeg) - 1._dl)
       else
         be_equil_calc = 0._dl
       end if


       return


       end function be_equil_calc


!-------------------------------------------------------------------     

       real(dl) function fd_equil_calc(eps,tratio,etadeg)


!------linkages.

!      called by - [subroutine] calc_weight_nuer1
!      calls     - none


!------remarks.

!      Calculates part of weight function for \nu-e scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: eps !input epsilon quantity
       real(dl), intent(in) :: tratio !input temperature ratio
       real(dl), intent(in) :: etadeg !input degeneracy parameter


!------procedure.


       if (eps.le.300._dl) then
         fd_equil_calc = 1._dl/(exp(eps*tratio - etadeg) + 1._dl)
       else
         fd_equil_calc = 0._dl
       end if


       return


       end function fd_equil_calc


!-------------------------------------------------------------------     

       real(dl) function gammar(arg)


!------linkages.

!      called by - [function] get_coul_corr
!      calls     - none


!------remarks.

!      Calculates log of real gamma function
!      arg > 0


       implicit none


!------throughput variables.

       real(dl), intent(in) :: arg !input argument


!------local variables.

       real(dl) part1, part2, sersum
       integer i


!------procedure.


       part1 = (arg + 0.5_dl)*log(arg + 5.5_dl) - (arg + 5.5_dl)
       sersum = 1._dl
       do i=1,6
         sersum = sersum + gamcoeffs(i)/(arg + real(i, kind=dl))
       end do
       part2 = log(stp*sersum/arg)

       gammar = part1 + part2
       

       return


       end function gammar


!-------------------------------------------------------------------     

       complex(dl) function gamcar(arg)


!------linkages.

!      called by - [function] get_coul_corr
!      calls     - none


!------remarks.

!      Calculates log of complex gamma function
!      Re(arg) > 0


       implicit none


!------throughput variables.

       complex(dl), intent(in) :: arg !input argument


!------local variables.

       complex(dl) part1, part2, sersum
       integer i


!------procedure.


       part1 = (arg + 0.5_dl)*log(arg + 5.5_dl) - (arg + 5.5_dl)
       sersum = complex(1._dl,0._dl)
       do i=1,6
         sersum = sersum + gamcoeffs(i)/(arg + real(i, kind=dl))
       end do
       part2 = log(stp*sersum/arg)

       gamcar = part1 + part2


       return


       end function gamcar


!-------------------------------------------------------------------     

       real(dl) function lep_to_xi(arg)


!------linkages.

!      called by - [program] main
!      calls     - none


!------remarks.

!      Calculates degeneracy parameter given lepton number
!      Assumes FD equilibrium


       implicit none


!------throughput variables.

       real(dl), intent(in) :: arg !input lepton number


!------local variables.

       real(dl) term1, y, y3


!------procedure.


       term1 = 12._dl*zeta3*arg
       if (arg.gt.0._dl) then
         y3 = term1/2._dl + 0.5_dl*sqrt(term1**2 + 4._dl*pi**2/27._dl)
         y = y3**(1.0/3.0)
         lep_to_xi = y - pi**2/3._dl/y
       else if (arg.lt.0._dl) then
         y3 = term1/2._dl - 0.5_dl*sqrt(term1**2 + 4._dl*pi**2/27._dl)
         y = -(abs(y3)**(1.0/3.0))
         lep_to_xi = y - pi**2/3._dl/y
       else
         lep_to_xi = 0._dl
       end if


       return


       end function lep_to_xi


       end module mainvar

