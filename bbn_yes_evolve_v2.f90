
       subroutine bbn_yes_evolve_v2(t,bbnvalues,bbnderivs,loop)

!------linkages.

!      called by - [subroutine] bbn_yes_rkck_v2
!      calls     - [subroutine] bbn_yes_therm_v2, trans_evolve,


!------remarks.

!      computes derivatives of
!       - temperature
!       - scale factor


!------modules.

       use bbnvar_v2
       use ratenp
       use nucsolver
#ifdef prllel
       use mpi
#endif       


       implicit none


       interface bbn_yes_therm_v2_interface
         subroutine bbn_yes_therm_v2(bbnvalues, tcmev, loop)
           use bbnvar_v2
           use bessel
           implicit none
           type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
           real(dl), intent(in) :: tcmev !comoving temperature in MeV
           integer, intent(in) :: loop
         end subroutine bbn_yes_therm_v2
       end interface bbn_yes_therm_v2_interface

       interface trans_evolve_interface
         subroutine trans_evolve(nutrans, doccfracdt, drhonudt, tcmev, tplmev, phie)
           use transvar
           use trans_nunu
           use trans_nunubar
           use trans_nuer1
           use trans_nuer2
           use trans_nuepma
#ifdef prllel
           use mpi
#endif
           implicit none
           type(nuvar), dimension(:,:), intent(in) :: nutrans !\nu variable
           real(dl), dimension(:,:,:), intent(out) :: doccfracdt !\nu occfrac derivatives
           real(dl), intent(out) :: drhonudt !change in \nu energy density
           real(dl), intent(in) :: tcmev !comoving temp. in MeV
           real(dl), intent(in) :: tplmev !Plasma temp. in MeV
           real(dl), intent(in) :: phie !e^- deg. parameter
         end subroutine trans_evolve
       end interface trans_evolve_interface

       !interface ratenuc_interface
       !  subroutine ratenuc(t9)
       !    use bbnvar_v2
       !    implicit none
       !    real(dl), intent(in) :: t9
       !  end subroutine ratenuc
       !end interface ratenuc_interface

       !save


!------throughput variables.

       real(dl), intent(in) :: t !time
       type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
       !type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
       type(bbnevolvar), intent(out) :: bbnderivs !bbn evolution variable derivatives
       integer, intent(in) :: loop


!------local variables.

       integer i,n,m !indices
       real(dl), allocatable, dimension(:,:,:) :: doccfracdt !\nu occfrac derivatives
       real(dl) sumy         !sum of abundances.
       real(dl) sumzy        !sum of charge*abundances.
       real(dl) sumdy        !sum of abundance flows.
       real(dl) summdy       !sum of (mass excess)*(abundance flows).
       real(dl) sumzdy       !sum of (charge)*(abundance flows).
       real(dl) dphdtpl       !d(phi e)/d(t9).
       real(dl) dphdln       !d(phi e)/d(h).
       real(dl) dphdzy       !d(phi e)/d(sumzy).
       real(dl) bar          !baryon density and pressure terms.
       real(dl) drhonudt
       real(dl) tpl, tplmev !plasma temp. in 10^9 K and MeV
       real(dl) hv, phie
       real(dl) asf !scale factor
       real(dl) tcm !comoving temperatures
       real(dl) :: tcmev !comoving temperatures
       real(dl) dtpl, dasfdt !derivatives
       real(dl) dhv                  !change in hv.
       real(dl) dphie                !change in chemical potential.
       real(dl) dsendt         !entropy per baryon time derivative
       real(dl) fori, revi
#ifdef prllel
       integer ierr
#endif


!------procedure.

       tpl = bbnvalues%tpl
       tplmev = boltzmann*tpl*1.e-06_dl
       hv = bbnvalues%hv
       phie = bbnvalues%phie
       asf = bbnvalues%asf

       !if (loop.eq.6) y0 = bbnvalues%y

       if (.not.allocated(doccfracdt)) then
         allocate(doccfracdt(size(bbnvalues%nuoccprob,1), &
                  2,size(bbnvalues%nuoccprob(1,1)%occfrac)))
       end if

!10----compute derivatives for abundances-------------

       rnb = 1._dl/asf**3

       if (tplmev.lt.8._dl) then
         tcm = (asfdec/asf)*tcmdec
       else
         tcm = tpl
         asfdec = asf
         tcmdec = tcm
       end if
       tcmev = boltzmann*tcm*1.e-06_dl

!------various thermodynamic quantities.

       call bbn_yes_therm_v2(bbnvalues, tcmev, loop)
       !15Apr2015 EG: Implement hubcalc & hubparm
       !hubcst = sqrt(8._dl/3._dl*pi*g*thm(10)) !expansion rate.
       hubcst = sqrt(8._dl/3._dl*pi/mpl**2*thm(10)/mev4tocgs)/hbar !expansion rate.
       !rhob = thm(9) !baryon mass density.

!------compute \nu energy transport.


       if (nurhoflag) then
         call trans_evolve(bbnvalues%nuoccprob,doccfracdt, &
              drhonudt, tcmev, tplmev, phie)
       else
         doccfracdt = 0._dl
         drhonudt = 0._dl
       end if

!------compute reaction rate coefficients.

       if (loop.eq.1) then
         !call ratenp_calc(tcmev,tcmev/tplmev,phie,cnorm,bbnvalues%nuoccprob(1,:),for1,rev1)
         call ratenp_calc(tcmev,tcmev/tplmev,phie,cnorm,bbnvalues%nuoccprob(1,:),for1,rev1 &
              ,ip,it)
         dyetest3dt = for1 - yetest3*(for1 + rev1)
         rhob1 = thm(9)
       else if (loop.eq.5) then
         !call ratenp_calc(tcmev,tcmev/tplmev,phie,cnorm,bbnvalues%nuoccprob(1,:),for2,rev2)
         call ratenp_calc(tcmev,tcmev/tplmev,phie,cnorm,bbnvalues%nuoccprob(1,:),for2,rev2 &
              ,ip,it)
         rhob2 = thm(9)
       end if
       !call ratenuc(tpl) !forward rate for reactions with a < 10.

!------solve coupled differential equations.

       !call sol_v2(bbnvalues%y,bbnderivs%y,dt,bbnvalues%tpl)
       !call sol_v2(bbnvalues%y,bbnderivs%y,dt,bbnvalues%tpl,loop)
       !dydt = bbnderivs%y
       !if (mbad.gt.0) return !abort in case matrix not invertible.

!------accumulate to get sum.

       sumy = 0._dl
       sumzy = 0._dl
       sumdy = 0._dl
       summdy = 0._dl
       sumzdy = 0._dl
       do i = 1,isize
         sumy = sumy + bbnvalues%y(i)           !sum of abundance.
         sumzy = sumzy + zm(i)*bbnvalues%y(i)     !sum of charge*abundance.
         sumdy = sumdy + dydt(i)        !sum of abundance flow.
         summdy = summdy + dm(i)*dydt(i) !sum of (mass excess)*(abundance flow).
         sumzdy = sumzdy + zm(i)*dydt(i) !sum of (charge)*(abundance flow).
       end do

!20----compute derivatives for temperature, and scale factor

       !dphdtpl = thm(12)*(-1.070e-04_dl*hv*sumzy/tpl - thm(11))
       dphdtpl = thm(12)*(-3._dl*convfact7*hv*sumzy/tpl - thm(11))
       !dphdln = -thm(12)*3.568e-05_dl*hv*sumzy
       dphdln = -thm(12)*convfact7*hv*sumzy
       !dphdzy = thm(12)*3.568e-05_dl*hv
       dphdzy = thm(12)*convfact7*hv
       !bar = 9.25e-05_dl*tpl*sumy + 1.388e-04_dl*tpl*sumdy/(3._dl*hubcst) &
       bar = 2._dl/3._dl*convfact6*tpl*sumy + convfact6*tpl*sumdy/(3._dl*hubcst) &
             + summdy/(3._dl*hubcst)
       !write (*,*) convfact6, convfact7

       dtpl = -(3._dl*hubcst)*(thm(1) + thm(3) + thm(4) + thm(7) + &
             thm(9)*bar + thm(6)*(dphdln + dphdzy*sumzdy/(3._dl*hubcst)) &
             + mev4tocgs*drhonudt/(3._dl*hubcst))/ &
             !(thm(2) + thm(5) + thm(6)*dphdtpl + thm(9)*1.388e-04_dl*sumy)   !(ref 1).
             (thm(2) + thm(5) + thm(6)*dphdtpl + thm(9)*convfact6*sumy)   !(ref 1).
       dlt9dt = dtpl/tpl

       dhv = -hv*((3._dl*hubcst) + 3._dl*dlt9dt)                    !(ref 2).

       dphie = dphdtpl*dtpl + dphdln*(3._dl*hubcst) + dphdzy*sumzdy  !(ref 3).

       dasfdt = asf*hubcst

       dsendt = -mev4tocgs*drhonudt/thm(9)*amu/tplmev !&
                !- convfact6*tpl*sumdy*amu/tplmev &
                !- summdy*amu/tplmev &
                !- thm(6)*(dphdln*3._dl*hubcst + dphdzy*sumzdy)/thm(9)*amu/tplmev

       if (loop.eq.1) then
         dspldt = dsendt
         dsnudt = 0._dl        
         do m=1,size(bbnvalues%nuoccprob, 1)
           do n=1,2
             do i=2,size(bbnvalues%nuoccprob(1,1)%occfrac)
               dsnudt = dsnudt &
                      + eps_bins%weights(i)*eps_bins%abscissas(i)**2 &
                      *doccfracdt(m,n,i) &
                      *log(1._dl/bbnvalues%nuoccprob(m,n)%occfrac(i) - 1._dl)
             end do !i
           end do !n
         end do !m
         dsnudt = mev4tocgs/2._dl/pi**2*tcmev**3*amu/thm(9)*dsnudt
         dstotdt = dspldt + dsnudt
       end if

!------msw quantities.
       if (loop.eq.1) then
         tcmmsw = tcmev
         yemsw = sum(zm*bbnvalues%y)
         rhoemsw = thm(4)/mev4tocgs
         nbarymsw = thm(9)/amu/mev4tocgs
         dtpldtmsw = dtpl*boltzmann*1.e-06_dl
         dyedtmsw = for1 - yemsw*(for1 + rev1)
         drhoedtmsw = thm(5)/mev3tocgs*dtpldtmsw  &
                 + thm(6)/mev4tocgs*dphie
         dnepdtmsw = tplmev**3/pi**2*2._dl*(thm(11)*dtpl + dphie/thm(12)) &
                    + 3._dl*nepdiffmsw*dtpl/tpl
         hubmsw = hubcst
       end if

       bbnderivs%tpl = dtpl
       bbnderivs%hv = dhv
       bbnderivs%phie = dphie
       !if (loop.eq.6) bbnvalues%y = y0
       !if (loop.eq.5) dydtstep = bbnderivs%y
       bbnderivs%y = 0._dl
       bbnderivs%asf = dasfdt
       bbnderivs%sen = dsendt
       do m=1,size(bbnderivs%nuoccprob, 1)
         do n=1,2
           bbnderivs%nuoccprob(m,n)%occfrac(:) = doccfracdt(m,n,:)
           bbnderivs%nuoccprob(m,n)%dilfact = 0._dl
           bbnderivs%nuoccprob(m,n)%mass = 0._dl
         end do !n
       end do !m

#ifdef prllel
         if (rank.eq.0) then
#endif
         !write (*,*) is, ip, it, loop, tplmev
         !write (*,*) tcmev, drhonudt, doccfracdt(2,1,:)
#ifdef prllel
         end if
#endif
       !if ((rank.eq.0).and.(it.ge.6).and.(ip.ge.133)) then
       !  write (101,*) ip, tplmev, tcmev, drhonudt
       !end if
       if (isnan(tcmev)) then
#ifdef prllel
         if (rank.eq.0) then
#endif
         write (*,*) is, ip, it, tplmev
         write (*,*) drhonudt, dtpl, dhv, dphie
#ifdef prllel
         end if
         call mpi_barrier(mpi_comm_world,ierr)
#endif
         stop
       end if

       return

!----------references-----------------------------------------------
!     1)  kawano, l., 1992, fermilab preprint fermilab-pub-92/04-a,
!         kellogg radiation lab preprint oap-714,
!         equation d.35.
!     2)  kawano, l., 1992, preprint, equation d.19.
!     3)  kawano, l., 1992, preprint, equation d.20.


       end subroutine bbn_yes_evolve_v2

