
       subroutine bbn_yes_launch_v2(t,dt,bbnvalues,obh2)


!------linkages.

!      called by - [subroutine] bbn_yes
!      calls     - [subroutine] cscgqf, nse_get, ratedecay


!------remarks.

!      Launches a bbn computation.


!------modules.

       use bbnvar_v2
       use bessel
       use nse
       use renorm
       use ratenp


       implicit none


!------throughput variables.

       real(dl), intent(out) :: t !time
       real(dl), intent(out) :: dt !time step
       type(bbnevolvar), intent(out) :: bbnvalues !bbn evolution variable
       real(dl) obh2 !baryon density


!------local variables.

       integer i,n,m        !indicies
       real(dl) hubcst0       !defined by hubcst=(8/3*pi*g*rho+lambda/3)^0.5
       real(dl) ggs           !initial statistical weight in entropy
       real(dl) sfsi
       real(dl) z             !defined by z = m(electron)*c**2/k*t9.
       character char1*200,char2*12,char3*12
       real(dl) ti !initial time
       real(dl) tcmev
       real(dl) enu
       real(dl) ee
       real(dl) phspint
       real(dl) coulcorr, zerorcorr
       real(dl) t9, tplmev
       real(dl) dmeeps2
       real(dl) dmedteps
       real(dl) dmgeps2
       real(dl) dmgdteps
       real(dl) x, dx, occprob, mt
       real(dl) glag1, glag3, glag4, glag7
       logical successflag


!------procedure.

!10----initialize flags and counters.


       is = 1 !first iteration coming up.
       ip = inc !set to maximum allowed # of iteration.
       it = 0 !no accumulation yet.
       mbad = 0 !no computational errors.


!20----settings.

!------computational settings.

       bbnvalues%tpl = t9i
       t9 = bbnvalues%tpl                   !initial temperature.
       tplmev = boltzmann*bbnvalues%tpl*1.e-06_dl
       !tcm =  bbnvalues%tpl !initial comoving temperature.
       tcmev = boltzmann*bbnvalues%tpl*1.e-06_dl
       !05Jun2015 EG: Need to insert calculation for time:
       ti = 1._dl/(const1*bbnvalues%tpl)**2 !initial time (ref 1).
       t = ti
       dt = dt1                    !initial time step.
       tcmdec = tcmev

!------model settings.

       !18Apr2015 EG: Need to change hubcst0 to something calculable
       hubcst0 = 608.5_dl     !initial expansion rate (s^-1) at t9~348?
       hubcst = hubcst0
       ggs = 11._dl/2._dl        !initial entropic statistical weight

       !Convert rho_{rad} into proper units at start of computation
       rhors = 232011._dl*rhors*tcmev**4*(11._dl/4._dl)**(4._dl/3._dl) !cgs


!------find constant in front of n<->p rates

       phspint = 0._dl
       sclgndra = nplgndr%abscissas
       sclgndrw = nplgndr%weights
       call cscgqf(1,sclgndra,sclgndrw,0._dl,deltamnp-xmelec,successflag)
       do i=1,size(sclgndra)
         enu = sclgndra(i)
         ee = deltamnp - enu
         coulcorr = get_coul_corr(-1._dl,sqrt(ee**2 - xmelec**2)/ee, &
                    sqrt(ee**2 - xmelec**2))
         zerorcorr = get_zero_rad_corr(sqrt(ee**2 - xmelec**2)/ee &
                                      , enu/xmelec, ee/xmelec)
         phspint = phspint + coulcorr*zerorcorr &
                   *enu**2*ee*sqrt(ee**2 - xmelec**2)*sclgndrw(i)
         !if (rank.eq.0) write (243, *) enu,ee,coulcorr, coulcorr*sqrt(ee**2 - xmelec**2), phspint
       end do
       cnorm=1._dl/mntau/phspint       !has units of 1/s/mev^5
       !if (rank.eq.0) write(*,*) cnorm


       !initialize reaction rates:
       f = 0._dl !forward rate coeff.
       r = 0._dl !reverse rate coeff.


!------Initialize flag for transport:

       nurhoflag = .false.

       !asf = 1._dl
       bbnvalues%asf = 1._dl
       asfdec = 1._dl


!30----compute initial abundances for neutron and proton

       !bbnvalues%y(1) = 1._dl/(exp(15.011_dl/t9+xi(1))+1._dl) !initial n abundance (ref 3).
       !bbnvalues%y(2) = 1._dl/(exp(-15.011_dl/t9-xi(1))+1._dl) !initial p abundance (ref 3).
       bbnvalues%y(1) = 1._dl/(exp(15.00849810551581_dl/t9+xi(1))+1._dl) !initial n abundance (ref 3).
       bbnvalues%y(2) = 1._dl/(exp(-15.00849810551581_dl/t9-xi(1))+1._dl) !initial p abundance (ref 3).

       yetest3 = 1._dl/(1._dl + bbnvalues%y(1)/bbnvalues%y(2))
       
!40----find ratio of baryon density to temperature cubed

       !21Aug2015 EG: Needs to be changed to incorporate QED
       !corrections:
       !z = 5.930_dl/t9            !inverse of temperature.
       !call get_eg_mass_renorm(tplmev,0._dl,glagbessel,dmeeps2,dmedteps,dmgeps2,dmgdteps)
       dmeeps2 = 0._dl
       dmedteps = 0._dl
       dmgeps2 = 0._dl
       dmgdteps = 0._dl
       !seninit = 121310._dl*ggs/2._dl/bbnvalues%hv
       glag1 = 0._dl
       glag3 = 0._dl
       glag4 = 0._dl
       glag7 = 0._dl
       mt = xmelec/tplmev
       do i=1,size(glagbessel%abscissas)
         x = glagbessel%abscissas(i)
         dx = glagbessel%weights(i)

         occprob = be_equil_calc(sqrt(x**2 + dmgeps2),1._dl,0._dl)

         glag1 = glag1 &
               + dx*x**2*sqrt(x**2 + dmgeps2)*occprob
         glag3 = glag3 &
               + dx*x**4/3._dl/sqrt(x**2 + dmgeps2)*occprob

         occprob = fd_equil_calc(sqrt(x**2 + mt**2 + dmeeps2),1._dl,0._dl)

         glag4 = glag4 &
               + dx*x**2*sqrt(x**2 + mt**2 + dmeeps2)*2._dl*occprob
         glag7 = glag7 &
               + dx*x**4/3._dl/sqrt(x**2 + mt**2 + dmeeps2)*2._dl*occprob
       end do
       glag1 = 1._dl/pi**2*glag1
       glag3 = 1._dl/pi**2*glag3
       glag4 = 1._dl/pi**2*glag4
       glag7 = 1._dl/pi**2*glag7
       thm(1) = mev4tocgs*tplmev**4*glag1
       thm(3) = mev4tocgs*tplmev**4*glag3
       thm(4) = mev4tocgs*tplmev**4*glag4
       thm(7) = mev4tocgs*tplmev**4*glag7
       !sfsi = (11/2)/2 = 11/4 with no QED corrections.
       !sfsi = (glag1 + glag3 + glag4 + glag7)/(4._dl*pi**2/45._dl) with
       !no transport corrections.
       sfsi = (glag1 + glag3 + glag4 + glag7)/(4._dl*pi**2/45._dl) !&
       if (transflag) sfsi = sfsi/1.0039771983175025_dl !Correction from trans when turning on at 8 MeV
       !if (transflag) sfsi = sfsi/1.0035742001227501_dl !Correction from trans with only annih
       !if (transflag) sfsi = sfsi/1.0007425811435222_dl !Correction from trans with only nue
       !bbnvalues%hv = 3.3683e+04_dl*eta1*2.75_dl !&
       !rhob0 = bbnvalues%hv*bbnvalues%tpl**3 !baryon density.
       rhob0 = mev4tocgs*convfact5*obh2*(bbnvalues%tpl/temp0k)**3*sfsi !baryon density.
       thm(9) = rhob0
       !bbnvalues%hv = 3.3683e+04_dl*eta1*sfsi !&
       bbnvalues%hv = rhob0/bbnvalues%tpl**3 !&
       !seninit = 121310._dl/bbnvalues%hv*sfsi
       seninit = (glag1+glag3+glag4+glag7)*tplmev**3*amu/(rhob0/mev4tocgs)
       !write (*,*) seninit, 121310._dl/bbnvalues%hv*11._dl/4._dl, seninit*(1._dl+25._dl*alphafs/22._dl/pi)
       bbnvalues%sen = seninit
       snu = 7._dl*pi**2/30._dl*tplmev**3*amu/(rhob0/mev4tocgs)
       stot = seninit + snu
       !chemical potential of electron (ref 5):
       z = sqrt(xmelec**2/tplmev**2 + dmeeps2)
       call bessel_eval(z)
       bbnvalues%phie = bbnvalues%hv*(1.784e-5_dl*bbnvalues%y(2))/(0.5_dl*z**3 &
               *(blz(1) - 2._dl*blz(2) + 3._dl*blz(3) - 4._dl*blz(4) + 5._dl*blz(5)))
       !write (*,*) rhob0, bbnvalues%sen, bbnvalues%hv, bbnvalues%phie
                 
 
!50--------set abundances for rest of nuclides-------------------

       !bbnvalues%y(3) = bbnvalues%y(1)*bbnvalues%y(2)*rhob0*exp(25.82_dl/t9)/(0.471e+10_dl*t9**1.5)!(ref 7).
       !do i = 4,isize
       !  bbnvalues%y(i) = ytmin !init abundances at beginning of iteration.
       !end do
       call nse_get(bbnvalues%tpl*boltzmann*1.e-06_dl,bbnvalues%y(1),bbnvalues%y(2),rhob0,bbnvalues%y)
       do i = 7,isize
         bbnvalues%y(i) = ytmin !init abundances at beginning of iteration.
       end do
       y0 = bbnvalues%y
       dydt = 0._dl

       heflag = .true.

       call ratedecay

       return


!----------references-----------------------------------------------
!     1) wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 44, equation a15.
!     2) coulomb correction obtained by dividing by correction factor fp(t9)
!        fp(t9) = 1 - 0.5(pi/(137<v>/c))
!        wagoner, r.v. 1973, ap. j. 179, page 358.
!     3) for the nondegenerate case:
!        wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 4, equation 3.
!        for the case with neutrino degeneracy:
!        beaudet,g. and goret,p., 1976, astron. & astrophys., 49,
!        page 417, equation 9.
!     4) wagoner, r.v. 1969, ap j. suppl. no. 162, 18, page 250, equation 4.
!        3.3683e+4 = mu(ng/t9**3) with mu the atomic mass, ng the
!        photon density.  2.75 is for the 11/4 factor difference
!        between the initial and final values of eta.
!     5) kawano, l., 1992, fermilab preprint fermilab-pub-92/04-a,
!        kellogg radiation lab preprint oap-714.
!        equation d.2.
!     6) wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 43, equation a4.
!        7.366 is used instead of 14.73 as the latter is the sum total
!        for 2 neutrino species.
!     7) initial deuterium abundance from nuclear statistical equilibrium
!        wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!        page 19, equation 17.


       end subroutine bbn_yes_launch_v2

