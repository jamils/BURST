
       module sdratiovar


!------linkages.

!      uses - [module] mainvar
!      used by - [subroutine] sdratio, sdratio_init, sdratio_launch
!                [module] xe_history, ang_diam_mod, opt_depth_mod, length_scale


!-----remarks.

!     Contains a list of variables and related functions


!-----modules.

       use mainvar


       implicit none


!------physical constants.

       real(dl) :: mpls = 1.220932_dl !Planck mass in 10^-22 MeV
       real(dl) :: mpro = 938.272046_dl !proton mass in MeV
       real(dl) :: melec = 0.510998928_dl !electron mass in MeV
       real(dl) :: ipoth1 = 13.60569253e-06_dl !H ionization potential in MeV
       real(dl) :: ipothe1 = 24.587401e-06_dl !HeI ionization potential in MeV
       real(dl) :: ipothe2 = 4._dl*13.60569253e-06_dl !HeII ionization potential in MeV
       real(dl) :: he1c2s = 20.6157734673e-06_dl !MeV
       real(dl) :: mh = 938.272046_dl + 0.510998928_dl - 13.60569253e-06_dl !H mass in MeV
       real(dl) :: spinh = 4._dl !internal spin partition function for H
       real(dl) :: spinp = 2._dl !internal spin partition function for free p^+
       real(dl) :: spine = 2._dl !internal spin partition function for free e^-
       real(dl) :: sigmat = 0.6652458734_dl/197.3269718_dl**2*1.e+02_dl !Thomson cross-section in units of MeV^-2:
       real(dl) :: alpha = 7.2973525698d-03 !fine-structure constant


!------input parameters.

       real(dl) reddcpl !Redshift at decoupling from Planck
       !real(dl) omegam !contribution from matter
       real(dl) obh2 !\Omega_b*h^2, baryons
       logical xeflag !flag to go past intnum for X_e evolution
       integer pastfact !factor to multiply intnum if xeflag = .true.
       logical neffevflag !flag to evolve N_{eff} as a function of scale factor


       real(dl) adcpl
       real(dl) masstol !tolerance to accept massl
       logical hierflag !Flag for neutrino mass hierarchy
       real(dl) summnulow !lowest value of \Sigma m_\nu
       real(dl) summnuhigh !highest value of \Sigma m_\nu
       real(dl), dimension(:), allocatable :: equilfrac !equilibrium occupation fractions.
       integer neffnum !array size
       real(dl) neffwidth
       real(dl), dimension(:,:), allocatable :: ratioarray
       real(dl), dimension(:), allocatable :: neffarray
       real(dl), dimension(:,:), allocatable :: ratioend !end point ratio values
       real(dl), dimension(:), allocatable :: ratiocalc !calculated ratio values for input cosmology
       real(dl), dimension(:), allocatable :: dratio !difference in ratio values
       real(dl), dimension(:), allocatable :: dneff1 !Derivative of neff wrt ratio at left endpoint
       real(dl), dimension(:), allocatable :: dneff2 !Deritative of neff wrt ratio at right endpoint
       integer neffevnum !dimension for number of points to store when neffevflag = .true.
       integer barind !index over barprec
       integer xieind !index over xieprec
       integer xixind !index over xixprec
       real(dl) neffrangel !lower limit on Neff range for theory calculations
       real(dl) neffrangeh !upper limit on Neff range for theory calculations

       integer barprec !number of grid points in obh2
       integer xieprec !number of grid points in xi_e
       integer xixprec !number of grid points in xi_x
       integer summnuprec !number of grid points in \Sigma m_\nu

       real(dl), dimension(:,:,:), allocatable :: runarray !stores abundance information for each run
#ifdef prllel
       real(dl), dimension(:,:,:), allocatable :: runarrayp !stores abundance information for each run
#endif
       logical writeflag


!------Neff numbers.

       real(dl) :: convfact1 = 197.3269718_dl*3.24077929e-38_dl !Conversion factor to go from MeV^-1 to Mpc
       !Conversion factor to use when using mpls:
       real(dl) :: convfact2 = 1.e-32_dl/(2.99792458e+03_dl)**2*(197.3269718_dl*3.24077929_dl)**2
       real(dl) neffval


!------Integration numbers.

       integer intnum !Number of integration points for integ_a
       integer rk4num !Number of integration points for rk4_stepper


!------Late universe numbers.

       real(dl) a0
       real(dl) nb0
       real(dl) ne0
       real(dl) rhob0
       real(dl) rhode
       real(dl) rhog0
       real(dl) temp0ev
       real(dl) t0mev


!------CMB numbers.

       !real(dl) adc


!------Chemical numbers.

       real(dl) ypxe !Y_p needed for the xeh module.
       real(dl), dimension(4) :: xeimax
       real(dl), dimension(3) :: xeiz = (/1._dl,1._dl,2._dl/) !No. of e^{-}'s for the xei


!------Hubble parameters.

       type(rhovar) :: hubparm !quantities needed for calculating the Hubble rate.


       contains


!--------------------------------------------------------------------

       real(dl) function hubcalcrho(rho) !This needs to be fixed to use changing rho


!------linkages.

!      called by - [subroutine] integ_a
!                  [function] xe_de_h, xe_de_he
!      calls     - none


!-----remarks.

!     Calculates the Hubble rate using the Friedmann equation.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: rho


!------procedure.

       hubcalcrho = sqrt(8._dl*pi/3._dl/mpl**2*rho)
     
       return
       

       end function hubcalcrho


!--------------------------------------------------------------------

       subroutine rk4_stepper(x,n,y,f,step,finegrid)
       

!------linkages.

!      called by - [function] xe_driver
!      calls     - [subroutine] f


!-----remarks.

!      Explicit fourth-order Runge-Kutta stepper algorithm.


       implicit none


!------throughput variables.

       real(dl) x
       integer, intent(in) :: n !Number of equations to integrate
       real(dl), dimension(n), intent(inout) :: y
       !real(dl), dimension(:), intent(inout) :: y
       real(dl), intent(in) :: step
       integer, intent(in) :: finegrid
       external :: f !subroutine for derivative(s)


!------local variables.

       real(dl) xin
       real(dl), dimension(n) :: yin
       real(dl), dimension(n) :: kout
       real(dl), dimension(n,4) :: k !array of derivatives
       real(dl) dx
       real(dl), dimension(n) :: dy

       integer j


!------procedure.

!------make fine grid.

       dx = step/real(finegrid, kind=dl)
       
!------Integrate through fine grid.

       do j=0,(rk4num-1)

         xin = x + real(j, kind=dl)*dx
         yin = y
         call f(xin,n,yin,kout)
         !call f(xin,yin,kout)
         k(:,1) = kout(:)

         xin = x + (0.5_dl + real(j, kind=dl))*dx
         yin(:) = y(:) + 0.5_dl*k(:,1)*dx
         call f(xin,n,yin,kout)
         !call f(xin,yin,kout)
         k(:,2) = kout(:)

         !xin for the third step is the same as the second step
         yin(:) = y(:) + 0.5_dl*k(:,2)*dx
         call f(xin,n,yin,kout)
         !call f(xin,yin,kout)
         k(:,3) = kout(:)

         xin = x + (1._dl + real(j, kind=dl))*dx
         yin(:) = y(:) + k(:,3)*dx
         call f(xin,n,yin,kout)
         !call f(xin,yin,kout)
         k(:,4) = kout(:)

         dy(:) =  dx/6._dl*(k(:,1) + 2._dl*k(:,2) + 2._dl*k(:,3) + k(:,4))

         y = y + dy

       end do
       
       
       return


       end subroutine rk4_stepper


!--------------------------------------------------------------------

       !subroutine rk6_stepper(x,n,y,f,step,yerr)
       subroutine rk6_stepper(x,y,f,step,yerr)
       

!------linkages.

!      called by - [function] xe_driver
!      calls     - [subroutine] f


!-----remarks.

!      Explicit sixth-order Runge-Kutta stepper algorithm.
!      Uses Cash-Karp method to take a step.

       implicit none

       interface recom_he_net_interface
         subroutine f(a,xe,dxeda)
           use mainvar
           implicit none
           real(dl), intent(in) :: a !scale factor
           !integer, intent(in) :: n !number of reactions in network
           !real(dl), dimension(n), intent(in) :: xe !free-electron fraction array
           real(dl), dimension(:), intent(in) :: xe !free-electron fraction array
           !real(dl), dimension(n), intent(out) :: dxeda !derivatives of xe
           real(dl), dimension(:), intent(out) :: dxeda !derivatives of xe
         end subroutine f
       end interface recom_he_net_interface




!------throughput variables.

       real(dl), intent(in) :: x
       !integer, intent(in) :: n !Number of equations to integrate
       !real(dl), dimension(n), intent(inout) :: y
       real(dl), dimension(:), intent(inout) :: y
       real(dl), intent(in) :: step
       !integer, intent(in) :: finegrid
       !real(dl), dimension(n), intent(inout) :: yerr
       real(dl), dimension(:), intent(inout) :: yerr
       external :: f !subroutine for derivative(s)


!------local variables.

       real(dl) xin
       !real(dl), dimension(n) :: yin
       real(dl), dimension(:), allocatable :: yin
       !real(dl), dimension(n) :: kout
       real(dl), dimension(:), allocatable :: kout
       !real(dl), dimension(n,6) :: k !array of derivatives
       real(dl), dimension(:,:), allocatable :: k !array of derivatives
       real(dl) dx

       integer j, l, m !indicies

       !real(dl), dimension(5) :: addx = (/0.2_dl, 0.3_dl, 0.6_dl, 1._dl, 0.875_dl/)
       real(dl), dimension(6) :: addx = (/0._dl, 0.2_dl, 0.3_dl, 0.6_dl, 1._dl, 0.875_dl/)
       real(dl), dimension(5,5) :: inty = transpose(reshape([ &
               0.2_dl,             0._dl,           0._dl,             0._dl,                0._dl, &
               3._dl/40._dl,       9._dl/40._dl,    0._dl,             0._dl,                0._dl, &
               0.3_dl,             -0.9_dl,         1.2_dl,            0._dl,                0._dl, &
               -11._dl/54._dl,     2.5_dl,          -70._dl/27._dl,    35._dl/27._dl,        0._dl, &
               1631._dl/55296._dl, 175._dl/512._dl, 575._dl/13824._dl, 44275._dl/110592._dl, 253._dl/4096._dl], &
               [5,5]))
       real(dl), dimension(6), parameter :: accumy = (/37._dl/378._dl, 0._dl, 250._dl/621._dl, &
                                                       125._dl/594._dl, 0._dl, 512._dl/1771._dl/)
       real(dl), dimension(6) :: diff54 = (/accumy(1) - 2825._dl/27648._dl, 0._dl, &
                                            accumy(3) - 18575._dl/48384._dl, accumy(4) - 13525._dl/55296._dl, &
                                            -277._dl/14336._dl, accumy(6) - 0.25_dl/)


!------procedure.

       if (.not.allocated(yin)) then
         allocate(yin(size(y)))
         allocate(kout(size(y)))
         allocate(k(size(y),6))
       end if


       dx = step


       do l=1,6
         xin = x + addx(l)*dx
         yin(:) = y(:)
         if (l.ge.2) then
           do m=1,(l-1)
             yin(:) = yin(:) + inty(l-1,m)*k(:,m)*dx
           end do
         end if
         call f(xin,yin,kout)
         k(:,l) = kout(:)
       end do

       do l=1,6
         y(:) = y(:) + accumy(l)*k(:,l)*dx
         yerr(:) = yerr(:) + diff54(l)*k(:,l)*dx
       end do
       
       
       return


       end subroutine rk6_stepper


!--------------------------------------------------------------------

       real(dl) function shcalc(a,neff)


!------linkages.

!      called by - [subroutine] recom, integ_a
!      calls     - none


!-----remarks.

!     Calculates the analytical expression for the sound horizon given N_{eff}.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), intent(in) :: neff !N_{eff}


!------local variables.

       real(dl) aeq !epoch of matter-radn. equality
       real(dl) keq !Wave number at matter-radn. equality
       real(dl) req !R at epoch of matter-radn. equality
       real(dl) rdc !R at epoch of photon decoupling
       real(dl) omh2 !(\Omega_b + \Omega_c)h^2


!------procedure.

       !Matter contribution:
       omh2 = obh2 + och2
       !epoch of matter-radiation equality in units of mpc:
       aeq = (1._dl + 7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval) &
             *8._dl*pi**3*t0mev**4/(mpls**2*45._dl*omh2)/convfact2
       !analytical soltion for sound horizon
       keq = sqrt(2._dl*omh2/aeq)/2.99792458e+03_dl
       req = 0.75_dl*obh2/omh2 &
             *(1._dl + 7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval)
       rdc = req*a/aeq
       shcalc = 2._dl/3._dl/keq*sqrt(6._dl/req)*log((sqrt(req + rdc) &
                + sqrt(1._dl + rdc))/(1._dl + sqrt(req)))


       return
       

       end function shcalc


       end module sdratiovar
