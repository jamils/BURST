
       module msw


!------remarks.

!      Module to implement Mikheyev-Smirnov-Wolfenstein effect


!------modules.

       use mainvar


       implicit none


       save


!------constants.

       real(dl) :: gfermi = 1.1663787e-11_dl !Fermi coupling constant (MeV^-2)
       real(dl) :: masswb = 8.0385e+04_dl !W boson mass (MeV)
       real(dl) :: masszb = 9.11876e+04_dl !Z boson mass (MeV)


!------sin and cosine of mixing angles.

       real(dl) s2theta, c2theta !\sin(2\theta), \cos(2\theta)


!------coefficients in quadratic eq. (all units MeV^2)

       real(dl) acoeff, bcoeff, ccoeff


!------epsres index.

       integer epsind !epsilon resonance index


!------quantities used by multiple subroutines.

       real(dl) etabg !baryon to photon ratio
       real(dl) rhonue !\rho_{\nue} + \rho_{\bnue} (MeV^4)
       real(dl) sumlepe !weighted sum of lepton numbers
       real(dl) deltam !mass difference (MeV)
       real(dl) depsdt !derivative of epsres (s^-1)

       real(dl), dimension(101) :: new_weights


       contains


!-----------------------------------------------------------

       subroutine msw_launch(nuparm)


!------linkages.

!      called by - [subroutine] bbn_yes
!      calls     - none


!------remarks.

!      Subroutine to launch MSW calculation.


       implicit none


!------throughput variables.

       type(nuvar), dimension(:,:), intent(in) :: nuparm


!------local variables.

       real(dl) s22t !\sin^2(2\theta)


!------procedure.

       !mixing angle in \sin^2(2\theta)
       if (size(nuparm, 1).eq.4) then
         !\sin(2\theta):
         s2theta = nuparm(4,1)%mix
         !\cos(2\theta):
         c2theta = sqrt(1._dl - s2theta**2)
         !epsilon resonance index:
         epsind = 0
         if (nuparm(4,1)%mass.ne.0._dl) then
           !mass difference (MeV)
           deltam = nuparm(4,1)%mass*1.e-06_dl
         else
           mswflag = .false.
         end if
       else !cannot proceed; will not call msw subroutines
         mswflag = .false.
       end if
 
       new_weights(:) = 1._dl
       new_weights(1) = 0.5_dl
       new_weights(101) = 0.5_dl


       return


       end subroutine msw_launch


!-----------------------------------------------------------

       subroutine msw_main(nuparm,tcm,tpl,ye,rhoe,nepdiff,nbary &
                          ,dtpldt,dyedt,drhoedt,dnepdt,hubrate,dtout)

!------linkages.

!      called by - [subroutine] bbn_yes_stepper_v2
!      calls     - [subroutine] msw_swap


!------remarks.

!      Master subroutine to calculate MSW quantities


       implicit none


!------throughput variables.

       type(nuvar), dimension(:,:), intent(inout) :: nuparm
       real(dl), intent(in) :: tcm !comoving temperature (MeV)
       real(dl), intent(in) :: tpl !plasma temperature (MeV)
       real(dl), intent(in) :: ye !electron fraction
       real(dl), intent(in) :: rhoe !e^- + e^+ energy density (MeV^4)
       real(dl), intent(in) :: nepdiff !n_{e^-} - n_{e^+} (MeV^3)
       real(dl), intent(in) :: nbary !number density of baryons (MeV^3)
       real(dl), intent(in) :: dtpldt !plasma temperature derivative (MeV/s)
       real(dl), intent(in) :: dyedt !ye derivative (s^-1)
       real(dl), intent(in) :: drhoedt !rhoe derivative (MeV^4/s)
       real(dl), intent(in) :: dnepdt !nepdiff derivative (MeV^3/s)
       real(dl), intent(in) :: hubrate !Hubble rate
       real(dl), intent(out) :: dtout !time step output


!------local variables.

       real(dl) epsres !epsilon resonance value
       real(dl) lzprob !Landau-Zeener probability
       real(dl) dtmdt !time derivative of mixing angle
                      !in instantaneous mass basis.
       real(dl) adbatparm !adiabaticity parameter
       logical successflag !flag which determines whether to proceed
                           !with msw subroutines

!------procedure.


       lzprob = 0._dl
       dtmdt = 0._dl
       adbatparm  = 0._dl

       call msw_epsres(tpl,nbary,tcm,ye,nuparm,rhoe,epsres,successflag)


       if (successflag) call &
          msw_epsind(epsres,epsind,size(nuparm(1,1)%occfrac),successflag)


       if (successflag) call msw_lzprob(epsres,hubrate,tpl,tcm,ye &
                  ,dtpldt,drhoedt,dyedt,lzprob,adbatparm,successflag)


       !if (successflag) call &
       !   msw_adiabaticity(epsres,tcm,nepdiff,dnepdt,hubrate,dtmdt,successflag)


       if (successflag) call msw_swap(nuparm,epsind,lzprob)


       if ((epsres.ne.0._dl).or.(depsdt.ne.0._dl)) then
         dtout = epsres/abs(depsdt)*0.1_dl
       else
         dtout = 1._dl/hubrate*1.e+06_dl
       end if



#ifdef prllel
       if (rank.eq.0) then
#endif
         write (102,*) tcm,epsres,lzprob,adbatparm
#ifdef prllel
       end if
#endif


       return


       end subroutine msw_main


!-----------------------------------------------------------

       subroutine msw_epsres(tpl,nbary,tcm,ye,nuparm,rhoe &
                            ,epsres,successflag)


!------linkages.

!      called by - [subroutine] msw_main
!      calls     - none


!------remarks.

!     Subroutine to find resonance epsilon


       implicit none


!------throughput variables.

       real(dl), intent(in) :: tpl !plasma temperature (MeV)
       real(dl), intent(in) :: nbary !number density of baryons (MeV^3)
       real(dl), intent(in) :: tcm !comoving temperature (MeV)
       real(dl), intent(in) :: ye !electron fraction
       type(nuvar), dimension(:,:), intent(in) :: nuparm
       real(dl), intent(in) :: rhoe !e^- + e^+ energy density (MeV^4)
       real(dl), intent(out) :: epsres !resonance epsilon
       logical, intent(out) :: successflag !flag indicating success


!------local variables.

       real(dl) ndensg !number density of photons (MeV^3)
       real(dl), dimension(3) :: ndensanu !active neutrino number density array (MeV^3)
       integer m !index


!------procedure.


       ndensg = 2._dl*zeta3/pi**2*tpl**3

       etabg = nbary/ndensg

       rhonue = 0._dl
       do m=1,3
         !ndensanu(m) = tcm**3/2._dl/pi**2*(sum(eps_bins%weights &
         ndensanu(m) = tcm**3/2._dl/pi**2*(sum(new_weights &
                                  *eps_bins%abscissas**2 &
                                  *nuparm(m,1)%occfrac) &
                              !-sum(eps_bins%weights &
                              -sum(new_weights &
                                  *eps_bins%abscissas**2 &
                                  *nuparm(m,2)%occfrac))
         if (m.eq.1) then
           !rhonue = tcm**4/2._dl/pi**2*(sum(eps_bins%weights &
           rhonue = tcm**4/2._dl/pi**2*(sum(new_weights &
                               *eps_bins%abscissas**3 &
                               *nuparm(m,1)%occfrac) &
                           !+sum(eps_bins%weights &
                           +sum(new_weights &
                               *eps_bins%abscissas**3 &
                               *nuparm(m,2)%occfrac))
         end if
       end do !m
        

       sumlepe = 2._dl*ndensanu(1) + ndensanu(2) + ndensanu(3)
       sumlepe = sumlepe/ndensg


!------coefficients for ax^2 + bx + c = 0

       acoeff = 2._dl*tcm**2*8._dl*sqrt(2._dl)*gfermi/3._dl &
                *(rhoe/masswb**2 + rhonue/masszb**2)
       bcoeff = -4._dl*sqrt(2._dl)*zeta3*gfermi*tcm*tpl**3/pi**2 &
                *(sumlepe + 1.5_dl*etabg*(ye - 1._dl/3._dl))
       ccoeff = c2theta*deltam**2


!------quadratic formula.

       if ((bcoeff**2/4._dl/acoeff/ccoeff).ge.1) then
         epsres = -bcoeff/2._dl/acoeff - 1._dl/2._dl/acoeff &
                  *sqrt(bcoeff**2 - 4._dl*acoeff*ccoeff)
         successflag = .true.
       else
         epsres = 0._dl
         successflag = .false.
       end if
       !if (rank.eq.0) write (102,*) tcm, sum(eps_bins%weights &
       !                           *eps_bins%abscissas**2 &
       !                           *nuparm(1,1)%occfrac), &
       !                       sum(eps_bins%weights &
       !                           *eps_bins%abscissas**2 &
       !                           *nuparm(1,2)%occfrac)

       if ((epsind.eq.6).or.(epsind.eq.5)) write(*,*) tcm, epsres, ndensg, ndensanu(:)

       return


       end subroutine msw_epsres


!-----------------------------------------------------------

       subroutine msw_epsind(epsres,epsind,occfracsize,successflag)


!------linkages.

!      called by - [subroutine] msw_main
!      calls     - none


!------remarks.

!      Finds epsilon index for swap


       implicit none


!------throughput variables.

       real(dl), intent(in) :: epsres !resonance epsilon
       integer, intent(inout) :: epsind !resonance bin number
       integer, intent(in) :: occfracsize !number of bins
       logical, intent(out) :: successflag !flag indicating success


!------local variables.

       integer testind


!------procedure.


       testind = int(epsres/maxepslinsave &
                 *real(occfracsize - 1, kind=dl) + 1._dl)

       if (testind.eq.(epsind + 1)) then
         epsind = testind
         successflag = .true.
       else
         successflag = .false.
       end if


       return


       end subroutine msw_epsind


!-----------------------------------------------------------

       subroutine msw_lzprob(epsres,hubrate,tpl,tcm,ye &
                  ,dtpldt,drhoedt,dyedt,lzprob,adbatparm,successflag)


!------linkages.

!      called by - [subroutine] msw_main
!      calls     - none


!------remarks.

!      Calculates Landau-Zeener probability


       implicit none


!------throughput variables.

       real(dl), intent(in) :: epsres !resonance epsilon
       real(dl), intent(in) :: hubrate !Hubble rate (s^-1)
       real(dl), intent(in) :: tpl !plasma temperature (MeV)
       real(dl), intent(in) :: tcm !comoving temperature (MeV)
       real(dl), intent(in) :: ye !electron fraction
       real(dl), intent(in) :: dtpldt !plasma temperature derivative (MeV/s)
       real(dl), intent(in) :: drhoedt !\rho_{e^-} + \rho_{e^+} temperature derivative (MeV^4/s)
       real(dl), intent(in) :: dyedt !ye derivative (s^-1)
       real(dl), intent(out) :: lzprob !probability of swap
       real(dl), intent(out) :: adbatparm !adiabaticity parameter
       logical, intent(out) :: successflag !flag to determine if lzprob calculated


!------local variables.

       real(dl) dsumlepedt !time derivative of sumlepe (s^-1)
       real(dl) detabgdt !time derivative of etabg (s^-1)
       real(dl) aderiv !derivative of acoeff (MeV^2/s)
       real(dl) bderiv !derivative of bcoeff (MeV^2/s)
       real(dl) vpot !potential (MeV)
       real(dl) dvpotdt !potential derivative (MeV/s)
       real(dl) scaleh !scale height (MeV^-1)
       real(dl) lresosc !resonance oscillation length (MeV^-1)


!------procedure.


       dsumlepedt = -3._dl*(hubrate + dtpldt/tpl)*sumlepe
       detabgdt = -3._dl*(hubrate + dtpldt/tpl)*etabg


!------derivatives of acoeff and bcoeff.

       aderiv = -2._dl*hubrate*acoeff &
               +2._dl*tcm**2*8._dl*sqrt(2._dl)*gfermi/3._dl &
                *(drhoedt/masswb**2 - 4._dl*hubrate*rhonue/masszb**2)
       bderiv = -4._dl*sqrt(2._dl)*zeta3*gfermi &
                 *(3._dl*tpl**2*dtpldt*tcm -tpl**3*tcm*hubrate)/pi**2 &
                 *(sumlepe + 1.5_dl*etabg*(ye - 1._dl/3._dl)) &
               - 4._dl*sqrt(2._dl)*zeta3*gfermi*tcm*tpl**3/pi**2 &
                *(dsumlepedt + 1.5_dl*(detabgdt*(ye - 1._dl/3._dl) &
                + etabg*dyedt))


!------Take time derivative of quadratic formula (ccoeff is constant).

       depsdt = -(acoeff*bderiv - bcoeff*aderiv)/2._dl/acoeff**2 &
               - (bcoeff*bderiv - 2._dl*ccoeff*aderiv)/2._dl/acoeff &
                /sqrt(bcoeff**2 - 4._dl*acoeff*ccoeff) &
               + sqrt(bcoeff**2 - 4._dl*acoeff*ccoeff)/2._dl/acoeff**2*aderiv


!------potential terms.

       vpot = -0.5_dl*(bcoeff/tcm + acoeff/tcm*epsres) !(MeV)
       dvpotdt = -0.5_dl*(bderiv/tcm + bcoeff/tcm*hubrate + aderiv/tcm*epsres &
                 +acoeff*(epsres/tcm*hubrate + depsdt/tcm)) !(MeV/s)


       !resonance oscillation length (MeV^-1):
       lresosc = 4._dl*pi*epsres*tcm/deltam**2/s2theta
       !scale height (MeV^-1):
       scaleh = vpot/abs(dvpotdt)*s2theta/c2theta/hbar
                            
!------adiabaticity parameter:

       if (lresosc.gt.0._dl) then
         adbatparm = 2._dl*pi*scaleh/lresosc
       else
         adbatparm = 0._dl
       end if

       if (adbatparm.gt.0._dl) then
         lzprob = exp(-pi/2._dl*adbatparm) !Landau-Zeener probability
         successflag = .true.
       else !issue in adbatparm
         lzprob = 0._dl
         successflag = .false.
       end if

       !if (rank.eq.0) write (102,*) tcm, scaleh, lresosc, lzprob

       return


       end subroutine msw_lzprob


!-----------------------------------------------------------

       subroutine msw_adiabaticity(epsres,tcm,nepdiff &
                                  ,dnepdt,hubrate,dtmdt,successflag)


!------linkages.

!      called by - [subroutine] msw_main
!      calls     - none


!------remarks.

!      Calculates adiabaticity parameter


       implicit none


!------throughput variables.

       real(dl), intent(in) :: epsres !resonance epsilon
       real(dl), intent(in) :: tcm !comoving temperature (MeV)
       real(dl), intent(in) :: nepdiff !n_{e^-} - n_{e^+} (MeV^3)
       real(dl), intent(in) :: dnepdt !nepdiff derivative (MeV^3/s)
       real(dl), intent(in) :: hubrate !Hubble rate
       real(dl), intent(out) :: dtmdt !time derivative of mixing angle
                                      !in instantaneous mass basis.
       logical, intent(out) :: successflag !flag indicating success


!------local variables.

       real(dl) edelta !\Delta quantity (MeV)
       real(dl) dedeltadt !derivative of edelta (MeV/s)
       real(dl) ampa !A amplitude (MeV)
       real(dl) dampadt !derivative of ampa (MeV/s)
       real(dl) thetam !mixing angle


!------procedure.


       edelta = deltam**2/2._dl/epsres/tcm
       dedeltadt = deltam**2/2._dl* &
                   (-depsdt/epsres**2/tcm + hubrate/epsres/tcm)

       ampa = sqrt(2._dl)*gfermi*nepdiff
       dampadt = sqrt(2._dl)*gfermi*dnepdt

       thetam = edelta*s2theta/(edelta*c2theta - ampa)

       dtmdt = 1._dl/(1._dl + thetam**2) &
               *((edelta*c2theta - ampa)*s2theta*dedeltadt &
               - edelta*s2theta*(c2theta*dedeltadt - dampadt)) &
               /(edelta*c2theta - ampa)**2

       if (dtmdt.gt.hubrate) then
         !no longer satisfies adiabaticity
         successflag = .false.
       else
         !adiabaticity still holds
         successflag = .true.
       end if

       !if (rank.eq.0) write (102,*) tcm, edelta, dedeltadt, ampa, dampadt, dtmdt

       return


       end subroutine msw_adiabaticity


!-----------------------------------------------------------

       subroutine msw_swap(nuparm,epsind,lzprob)


!------linkages.

!      called by - [subroutine] msw_main
!      calls     - none


!------remarks.

!      Swaps occupation probability.


       implicit none


!------throughput variables.

       type(nuvar), dimension(:,:), intent(inout) :: nuparm !neutrino variable
       integer, intent(in) :: epsind !resonance bin number
       real(dl), intent(in) :: lzprob !LZ probability


!------local variables.

       logical testind


!------procedure.


       !swap:
       nuparm(4,1)%occfrac(epsind) = nuparm(1,1)%occfrac(epsind) &
                                     *(1._dl - lzprob)
       nuparm(1,1)%occfrac(epsind) = nuparm(1,1)%occfrac(epsind) &
                                     *lzprob


       return


       end subroutine msw_swap


       end module msw
