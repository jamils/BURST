
       module bbnvar_v2


!------linkages.

!      uses - [module] mainvar
!      used by - [subroutine] bbn, accum, derivs, driver, nse, ratedecay, ratenp, start, therm, ratenuc
!                [module] bessel, nucsolver


!-----remarks.

!     Contains a list of variables and related functions


!-----modules.

       use mainvar


       implicit none


!------numbers from Kawano.

       integer, parameter :: ir=1             !input unit number. 'r'=read
       integer, parameter :: iw=1             !output unit number. 'w'=write
       integer, parameter :: nrec=34     !number of nuclear reactions.
       integer, parameter :: nnuc=9    !number of nuclides in calculation.
!       integer, parameter :: itmax=400 !maximum number of entries in output arrays.
       integer, parameter :: itmax=40000 !maximum number of entries in output arrays.
       real(dl) dtlow !lower limit on size of time step.
       real(dl), parameter :: const1=0.09615  !relation between time and temp.        
       real(dl), parameter :: g=6.6700e-8 !gravitational constant.
       real(dl), parameter :: q=2.531  !(mass(neutron)-mass(proton))/m(electron)
       integer, parameter :: iter=50          !number of gaussian quads.
       integer, parameter :: mord=1    !higher order in correction.
       real(dl) epstol !tolerance for convergence (.ge. 1.e-7).
       real(dl), parameter :: epschg = 100._dl
       real(dl), parameter :: epsmax = 100._dl
       real(dl), parameter :: convfact3 = 3._dl/8._dl/pi/amu*mpl**2/(2997.92458_dl)**2 &
                                          *(197.3269718_dl*3.24077929e-38_dl)**2 &
                                          *pi**2/zeta3/2._dl*1.e+18_dl
       real(dl), parameter :: convfact4 = 32._dl/135._dl*pi**3*amu/mpl**2*(2997.92458_dl)**2 &
                                          /(197.3269718_dl*3.24077929e-38_dl)**2
       real(dl), parameter :: convfact5 = 3._dl/8._dl/pi*mpl**2/(2997.92458_dl)**2 &
                                          *(197.3269718_dl*3.24077929e-38_dl)**2
       real(dl), parameter :: convfact6 = 3._dl/2._dl*boltzmann/amu*1.e-06_dl
       real(dl), parameter :: convfact7 = pi**2/2._dl*(xmelec/(boltzmann*1.e-06))**3*6.022140857e+23_dl*197.3269788**3*1.e-39
       real(dl) :: mev4tocgs =2.320115553791773e+05_dl
       real(dl) :: mev3tocgs = 1.999320693343379e+04_dl



!------computational parameters.

       real(dl) cy          !time step limiting constant on abundances.
       real(dl) ct          !time step limiting constant on temperature.
       real(dl) t9i         !initial temperature (in 10**9 k).
       real(dl) t9f         !final temperature (in 10**9 k).
       real(dl) ytmin       !smallest abundances allowed.
       integer inc         !accumulation increment.
       real(dl) phiei


!------counters and flags.

       integer is            !# total time steps for particular run.
       integer ip            !# time steps after outputting a line.
       integer it            !# times accumulated in output buffer.
       integer mbad          !indicates if gaussian elimination fails.


!------early universe model parameters.

       real(dl) mntau !neutron lifetime


!------output arrays.

       real(dl), dimension(itmax,nnuc) :: xout     !nuclide mass fractions.
       real(dl), dimension(itmax,8) :: thmout      !thermodynamic variables.
       real(dl), dimension(itmax) :: t9out        !temperature (in units of 10**9 k).
       real(dl), dimension(itmax) :: tout          !time.
       real(dl), dimension(itmax) :: dtout         !time step.
       real(dl), dimension(itmax) :: senout        !entropy ber baryon.
       real(dl), dimension(itmax) :: hubout        !expansion rate.
       real(dl), dimension(itmax) :: tcmevout      !co-moving temperature


!------reaction rates.

       real(dl), dimension(nrec) :: f              !forward reaction rate coefficients.        
       real(dl), dimension(nrec) :: r              !reverse reaction rate coefficients.        


!------reaction parameters.
       integer, dimension(nrec) :: iform          !reaction type code (1-11).
       integer, dimension(nrec) :: ii             !incoming nuclide type (1-28).
       integer, dimension(nrec) :: jj             !incoming light nuclide type (1-6).        
       integer, dimension(nrec) :: kk             !outgoing light nuclide type (1-6).
       integer, dimension(nrec) :: ll             !outgoing nuclide type (1-28).
       real(dl), dimension(nrec) :: rev            !reverse reaction coefficient.
       real(dl), dimension(nrec) :: q9             !energy released in reaction.


!------run option.

       integer irun                 !run network size.
       integer isize                !number of nuclides in computation.
       integer jsize                !number of reactions in computation.


!------variational parameters.

       real(dl) dt1                  !initial time step.
       real(dl) eta1                 !baryon-to-photon ratio.


!------entropy generation quantities.

       real(dl) seninit   !initial entropy per baryon
       real(dl) sen       !entropy per baryon


!------original values of entropy.

       real(dl) sen0          !value of entropy per baryon
       real(dl) dsendt0       !derivative value of entropy per baryon


!------evolution parameters.

       !real(dl) t9          !temperature (in units of 10**9 k).
       !real(dl) hv          !defined by hv = m(atomic)n(baryon)/t9**3.
       !real(dl) phie        !chemical potential for electron.
       !real(dl), dimension(nnuc) :: y     !relative number abundances.


!------evolution parameters (derivatives).

       !real(dl) dt9                  !change in temperature.


!------time and time step variables.

       !real(dl) ti                  !initial time
       !real(dl) t                    !time.
       !real(dl) dt             !time step, limit on time step.
       real(dl) dlt9dt               !(1/t9)*d(t9)/d(t).
       real(dl), dimension(nnuc) :: dydtstep
       real(dl), dimension(nnuc) :: dydt
       real(dl), dimension(nnuc) :: y0
       real(dl) for1,rev1,for2,rev2


!------neutrino parameters.
       !real(dl) t9mev !plasma temperature in mev.
       !real(dl) tcm !co-moving temperature.
       !real(dl) tcmev !co-moving temperature in mev.
       real(dl) tmevi !initial temperature in mev.
       real(dl) cnorm !normalization constant for n<->p rates
       !real(dl) rhonu !neutrino energy density.
       logical nurhoflag !if .true., numerically integrates energy densities.
       !real(dl), dimension(:), allocatable :: eqocc !FD equil \nu occ. probs.


!------energy densities.

       real(dl) rhob0 !initial baryon mass density.
       !real(dl) rhob !baryon mass density.
       real(dl) rhob1 !baryon mass density.
       real(dl) rhob2 !baryon mass density.
       real(dl) rnb !baryon mass density (ratio to init value).
       !real(dl) asf !scale factor variable before e^\pm annihilation
       !real(dl) dasfdt !time derivative of scale factor variable before e^\pm annihilation
       real(dl) rhors !initial radiation energy density.


!------dynamic variables.

       real(dl), dimension(14) :: thm              !thermodynamic variables.
       real(dl) hubcst               !expansion rate.

       
!------flag for computing BBN.

       logical bbncompflag


!------Gauss-Laguerre weights and abscissas for bbn_therm.

       type(bin_scheme) glagbessel


!------BBN evolution variable.

       type bbnevolvar

         real(dl) :: tpl
         real(dl) :: hv
         real(dl) :: phie
         real(dl), dimension(nnuc) :: y
         real(dl) :: asf
         real(dl) :: sen
         type(nuvar), dimension(:,:), pointer :: nuoccprob

       end type bbnevolvar


!------misc. parameters.

       integer imod,jmod      !indicies for module


!------testing parameters.

       real(dl) yetest3
       real(dl) dyetest3dt
       real(dl) ndensem, ndensep
       real(dl) dspldt, dsnudt, dstotdt
       real(dl) snu, stot
       logical heflag


!------msw quantities.

       real(dl)  tcmmsw !comoving temperature (MeV)
       real(dl)  yemsw !electron fraction
       real(dl)  rhoemsw !\rho_{e^-} + \rho_{e^+} (MeV^4)
       real(dl)  nbarymsw !n_b (MeV^3))
       real(dl)  nepdiffmsw !n_{e^-} - n_{e^+} (MeV^3)
       real(dl)  dtpldtmsw !derivative of plasma temp. (MeV/s)
       real(dl)  dyedtmsw !derivative of ye (s^-1)
       real(dl)  drhoedtmsw !derivative of rhoe (MeV^4/s)
       real(dl)  dnepdtmsw !derivative of nepdiff (MeV^3/s)
       real(dl)  hubmsw !Hubble rate (s^-1)


!------data

!    nuclide and corresponding number
!    --------------------------------
!    1) n         7) li6
!    2) p         8) li7
!    3) h2        9) be7
!    4) h3
!    5) he3
!    6) he4

!------nuclide data.

       real(dl), dimension(nnuc) :: am = & !atomic mass number of nuclide.
                                       (/1._dl,1._dl,2._dl,3._dl,3._dl,4._dl,6._dl,7._dl,7._dl/)

       real(dl), dimension(nnuc) :: zm = & !charge of nuclide (atomic number).
                                       (/0._dl,1._dl,1._dl,1._dl,2._dl,2._dl,3._dl,3._dl,4._dl/)

       real(dl), dimension(nnuc) :: dm = & !mass excess of nuclide.
                                       (/0.008665_dl,0.007825_dl,0.014102_dl,0.016050_dl, &
                                       0.016030_dl,0.002603_dl,0.015125_dl,0.016004_dl,0.016929_dl/)

     
!------reaction parameters values.

       real(dl), dimension(nrec,8) :: reacpr       !reaction parameters.

!------reaction rate coefficients (ref 1).

       data ((reacpr(imod,jmod),jmod=1,8),imod=1,11) / &
!              reac# type n1 n2 n3 n4 rev-coeff q-value
!              ----  ---- -- -- -- -- --------- -------
                   1.,1., 1.,0.,0., 2., 0.0  ,   0.0 , &    !n->p
                   2.,1., 4.,0.,0., 5., 0.0  ,   0.0 , &    !h3->he3
                   3.,4.,10.,0.,0., 6., 0.0  ,   0.0 , &    !li8->2he4
                   4.,1.,16.,0.,0.,17., 0.0  ,   0.0 , &    !b12->c12
                   5.,1.,21.,0.,0.,22., 0.0  ,   0.0 , &    !c14->n14
                   6.,4.,11.,0.,0., 6., 0.0  ,   0.0 , &    !b8->2he4
                   7.,1.,15.,0.,0.,14., 0.0  ,   0.0 , &    !c11->b11
                   8.,1.,18.,0.,0.,17., 0.0  ,   0.0 , &    !n12->c12
                   9.,1.,20.,0.,0.,19., 0.0  ,   0.0 , &    !n13->c13
                  10.,1.,23.,0.,0.,22., 0.0  ,   0.0 , &    !o14->n14
                  11.,1.,25.,0.,0.,24., 0.0  ,   0.0 /      !o15->n15
       data ((reacpr(imod,jmod),jmod=1,8),imod=12,22) / &
!              reac# type n1 n2 n3 n4 rev-coeff q-value
!              ----  ---- -- -- -- -- --------- -------
                  12.,2., 2.,1.,0., 3., 0.471,  25.82, &    !h(n,g)h2
                  13.,2., 3.,1.,0., 4., 1.63 ,  72.62, &    !h2(n,g)h3
                  14.,2., 5.,1.,0., 6., 2.61 , 238.81, &    !he3(n,g)he4         
                  15.,2., 7.,1.,0., 8., 1.19 ,  84.17, &    !li6(n,g)li7
                  16.,3., 5.,1.,2., 4., 1.002,   8.863,&    !he3(n,p)h3
                  17.,3., 9.,1.,2., 8., 0.998,  19.081,&    !be7(n,p)li7
                  18.,3., 7.,1.,4., 6., 1.070,  55.494,&    !li6(n,a)h3
                  19.,5., 9.,1.,0., 6., 4.70 , 220.39, &    !be7(n,a)he4
                  20.,2., 3.,2.,0., 5., 1.63 ,  63.750,&    !h2(p,g)he3
                  21.,2., 4.,2.,0., 6., 2.61 , 229.932,&    !h3(p,g)he4
                  22.,2., 7.,2.,0., 9., 1.19 ,  65.054/     !li6(p,g)be7
       data ((reacpr(imod,jmod),jmod=1,8),imod=23,34) / &
!              reac# type n1 n2 n3 n4 rev-coeff q-value
!              ----  ---- -- -- -- -- --------- -------
                  23.,3., 7.,2.,5., 6., 1.07 ,  46.631, &   !li6(p,a)he3
                  24.,5., 8.,2.,0., 6., 4.69 , 201.291, &   !li7(p,a)he4
                  25.,2., 6.,3.,0., 7., 1.53 ,  17.118, &   !h2(a,g)li6
                  26.,2., 6.,4.,0., 8., 1.11 ,  28.640, &   !h3(a,g)li7
                  27.,2., 6.,5.,0., 9., 1.11 ,  18.423, &   !he3(a,g)be7
                  28.,6., 3.,0.,1., 5., 1.73 ,  37.935, &   !h2(d,n)he3
                  29.,6., 3.,0.,2., 4., 1.73 ,  46.798, &   !h2(d,p)h3
                  30.,3., 4.,3.,1., 6., 5.54 , 204.117, &   !h3(d,n)he4
                  31.,3., 5.,3.,2., 6., 5.55 , 212.980, &   !he3(d,p)he4
                  32.,11.,5.,0.,2., 6., 3.39 , 149.230, &   !he3(he3,2p)he4      
                  33.,9., 8.,3.,1., 6., 9.95 , 175.476, &   !li7(d,na)he4
                  34.,9., 9.,3.,2., 6., 9.97 , 194.557/     !be7(d,pa)he4


!------mass values.

       real(dl) :: mpro = 938.272046_dl !proton mass in MeV
       real(dl) :: mneu = 938.272046_dl + 1.29333217_dl !neutron mass in mev


!------spin values.

       real(dl) :: spinp = 2._dl !proton spin degeneracy
       real(dl) :: spinn = 2._dl !neutron spin degeneracy


!------constants

!------comoving temperature quantities

       real(dl) asfdec
       real(dl) tcmdec

       contains


!-----------------------------------------------------------------------

       real(dl) function xintd(xlow,xhi,func,nq,opt1)


!------remarks.

!      Computes integral of func using 6-point Gauss-Legendre method.


       implicit none


!-----throughput variables.

       real(dl), intent(in) :: xlow !lower limit.
       real(dl), intent(in) :: xhi !upper limit.
       integer, intent(in) :: nq !number of six point gaussian quads.
       real(dl), optional, intent(in) :: opt1 !optional value


!------external functions.

       real(dl), external :: func !function to integrate


!------local variables.

       real(dl)	x !integration point
       real(dl) dist !size of quad interval.
       real(dl) cent !center of quad interval.
       real(dl) sumterm !summation of terms.
       integer nint !interval number.
       integer npnt !point number.


!----------abscissas and weight factors.

       integer, parameter :: np = 6 !6 point gaussian integration.       

       real(dl), dimension(np) :: u = (/-.93246951420315,-.66120938646627,-.23861918608320, &
               .23861918608320, .66120938646627, .93246951420315/)

       real(dl), dimension(np) :: w = (/.17132449237917,.36076157304814,.46791393457269, &
               .46791393457269,.36076157304814,.17132449237917/)


!------procedure.

!10----do integration.
       
       
       sumterm   = 0._dl

       dist  = (xhi - xlow)/real(nq, kind=dl) !size of quad interval.
       do nint = 1,nq
         cent = xlow + (real(nint, kind=dl) - 0.5_dl)*dist  !center of interval.
         do npnt = 1,np
           x = cent + 0.5_dl*dist*u(npnt) !integration point.
           if (present(opt1)) then
             sumterm = sumterm + func(x,opt1)*w(npnt)      !add up sum.
           else
             sumterm = sumterm + func(x)*w(npnt)      !add up sum.
           end if
         end do
       end do

!20----get integral value.

       xintd = sumterm*dist*0.5_dl !do integral.


       end function xintd


!-----------------------------------------------------------------------

       real(dl) function funcnu(eps,deg)


!------remarks.

!      Function for energy density of neutrino


       implicit none


!-----throughput variables.

       real(dl), intent(in) :: eps !epsilon value
       real(dl), intent(in) :: deg !nu degeneracy value


!------procedure.

       funcnu = eps**3/(1._dl + exp(eps - deg))


       end function funcnu


!-----------------------------------------------------------------------

       real(dl) function squarert(arg,tol)


!------remarks.

!      Alternative square root function to avoid imaginary numbers.


       implicit none


!-----throughput variables.

       real(dl), intent(in) :: arg         !argument and tolerance for square root function
       real(dl), intent(in) :: tol          !argument and tolerance for square root function

!------procedure.       
       
       if (abs(arg).lt.tol) then
         squarert = 0.0
       else
         squarert = arg**0.5         
       end if


       end function squarert


!-----------------------------------------------------------------------

       subroutine bbn_no_addtype_v2(bbn1,bbn2,realnum)
       

!------linkages.

!      called by - [subroutine] bbn_no_rkck_v2
!      calls     - none


!------remarks.

!      Subroutine to add type(bbnevolvar) to similar type


       !use bbnvar_v2


       implicit none


!------throughput variables.

       type(bbnevolvar), intent(inout) :: bbn1 !first bbn evolution variable
       type(bbnevolvar), intent(in) :: bbn2 !second bbn evolution variable
       real(dl), intent(in) :: realnum !number to multiply bbn2


!------local variables.

       integer m,n !indices



!------procedure.


       bbn1%tpl = bbn1%tpl + realnum*bbn2%tpl
       bbn1%asf = bbn1%asf + realnum*bbn2%asf

       !if (is.eq.1) write (*,*) bbn2%nuoccprob(1,1)%occfrac(1:3)

       do m=1,size(bbn1%nuoccprob,1)
         do n=1,2
           bbn1%nuoccprob(m,n)%occfrac(:) = bbn1%nuoccprob(m,n)%occfrac(:) &
                                       + realnum*bbn2%nuoccprob(m,n)%occfrac(:) 
         end do
       end do !m

       !if (is.eq.1) write (*,*) bbn2%nuoccprob(1,1)%occfrac(1:3)


       return


       end subroutine bbn_no_addtype_v2


!-----------------------------------------------------------------------

       subroutine bbn_yes_addtype_v2(bbn1,bbn2,realnum)
       

!------linkages.

!      called by - [subroutine] bbn_yes_rkck_v2
!      calls     - none


!------remarks.

!      Subroutine to add type(bbnevolvar) to similar type


       implicit none


!------throughput variables.

       type(bbnevolvar), intent(inout) :: bbn1 !first bbn evolution variable
       type(bbnevolvar), intent(in) :: bbn2 !second bbn evolution variable
       real(dl), intent(in) :: realnum !number to multiply bbn2


!------local variables.

       integer m,n !indices



!------procedure.


       bbn1%tpl = bbn1%tpl + realnum*bbn2%tpl
       bbn1%hv = bbn1%hv + realnum*bbn2%hv
       bbn1%phie = bbn1%phie + realnum*bbn2%phie
       bbn1%asf = bbn1%asf + realnum*bbn2%asf
       bbn1%sen = bbn1%sen + realnum*bbn2%sen

       !if (is.eq.1) write (*,*) bbn2%nuoccprob(1,1)%occfrac(1:3)

       do m=1,size(bbn1%nuoccprob,1)
         do n=1,2
           bbn1%nuoccprob(m,n)%occfrac(:) = bbn1%nuoccprob(m,n)%occfrac(:) &
                                       + realnum*bbn2%nuoccprob(m,n)%occfrac(:) 
         end do
       end do !m

       !if (is.eq.1) write (*,*) bbn2%nuoccprob(1,1)%occfrac(1:3)


       return


       end subroutine bbn_yes_addtype_v2


!-----------------------------------------------------------------------

       subroutine bbn_no_assigntype_v2(bbn1,bbn2)
       

!------linkages.

!      called by - [subroutine] bbn_no_rkck_v2
!      calls     - none


!------remarks.

!      Subroutine to assign bbn2 values to bbn1


       !use bbnvar_v2


       implicit none


!------throughput variables.

       type(bbnevolvar), intent(inout) :: bbn1 !first bbn evolution variable
       type(bbnevolvar), intent(in) :: bbn2 !second bbn evolution variable


!------local variables.

       integer m,n !indices



!------procedure.


       bbn1%tpl = bbn2%tpl
       bbn1%asf = bbn2%asf


       do m=1,size(bbn1%nuoccprob,1)
         do n=1,2
           bbn1%nuoccprob(m,n)%dilfact = bbn2%nuoccprob(m,n)%dilfact
           bbn1%nuoccprob(m,n)%mass = bbn2%nuoccprob(m,n)%mass
           bbn1%nuoccprob(m,n)%occfrac(:) = bbn2%nuoccprob(m,n)%occfrac(:)
         end do
       end do !m


       return


       end subroutine bbn_no_assigntype_v2


!-----------------------------------------------------------------------

       subroutine bbn_yes_assigntype_v2(bbn1,bbn2)
       

!------linkages.

!      called by - [subroutine] bbn_yes_rkck_v2
!      calls     - none


!------remarks.

!      Subroutine to assign bbn2 values to bbn1


       implicit none


!------throughput variables.

       type(bbnevolvar), intent(inout) :: bbn1 !first bbn evolution variable
       type(bbnevolvar), intent(in) :: bbn2 !second bbn evolution variable


!------local variables.

       integer m,n !indices



!------procedure.


       bbn1%tpl = bbn2%tpl
       bbn1%hv = bbn2%hv
       bbn1%phie = bbn2%phie
       bbn1%y = bbn2%y
       bbn1%asf = bbn2%asf
       bbn1%sen = bbn2%sen


       do m=1,size(bbn1%nuoccprob,1)
         do n=1,2
           bbn1%nuoccprob(m,n)%dilfact = bbn2%nuoccprob(m,n)%dilfact
           bbn1%nuoccprob(m,n)%mass = bbn2%nuoccprob(m,n)%mass
           bbn1%nuoccprob(m,n)%occfrac(:) = bbn2%nuoccprob(m,n)%occfrac(:)
         end do
       end do !m


       return


       end subroutine bbn_yes_assigntype_v2


       end module bbnvar_v2
