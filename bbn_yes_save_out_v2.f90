
       subroutine bbn_yes_save_out_v2(t,bbnvalues,dt,eqocc)


!------linkages.

!      called by - [subroutine] bbn_yes_driver_v2
!      calls     - [subroutine] nse_det, ye_det


!------remarks.

!      output accumulator.


!------modules.

       use bbnvar_v2
       use nse
       use ye
#ifdef prllel
       use mpi
#endif


       implicit none


       interface bbn_yes_save_state_v2_interface
         subroutine bbn_yes_save_state_v2(t,bbnvalues,dt)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: t !time
           type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
           real(dl), intent(in) :: dt !time step
         end subroutine bbn_yes_save_state_v2
       end interface bbn_yes_save_state_v2_interface


!------throughput variables.

       real(dl), intent(in) :: t !time
       type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
       real(dl), intent(in) :: dt !time step
       real(dl), dimension(:), intent(in) :: eqocc !equilibrium occupation probabilities


!------local variables.

       integer i,j !indicies
       real(dl) tcmev, tplmev !temperatures
       real(dl) tcm
       real(dl) alln, allp, allnp
       character char1*200,char2*12
       real(dl) sentest,sentest2
       real(dl) barnum
#ifdef prllel
       integer ierr
#endif


!------procedure.

       it = it + 1                  !set up accumulation counter.


!40----set up output variables.

!------divide number fraction by that of proton.
       do i = 1,isize
         xout(it,i) = bbnvalues%y(i)/bbnvalues%y(2)
       end do
       xout(it,2) = bbnvalues%y(2)*am(2)      !exception for proton.
       xout(it,6) = bbnvalues%y(6)*am(6)      !exception for helium.

       tplmev = boltzmann*bbnvalues%tpl*1.e-06_dl
       !08Jun2015 EG: approximate for now:
       if (tplmev.lt.8._dl) then
         tcm = (asfdec/bbnvalues%asf)*tcmdec
       else
         tcm = bbnvalues%tpl
       end if
       tcmev = boltzmann*tcm*1.e-06_dl

!------calculate n/p ratio for all n and p
       alln = 0._dl
       allp = 0._dl
       do i=1,isize
         alln = alln + (am(i) - zm(i))*bbnvalues%y(i)
         allp = allp + zm(i)*bbnvalues%y(i)
       end do
       allnp = alln/allp

!------relabel temperature, time, thermodynamic variables, etc.
       t9out(it) = bbnvalues%tpl !temperature.
       tcmevout(it) = tcmev            !co-moving temperature.
       tout(it) = t             !time.
       thmout(it,1) = thm(1)       !rho photon.
       thmout(it,2) = thm(4)        !rho electron.
       thmout(it,3) = thm(8)        !rho neutrino.
       thmout(it,4) = thm(9)        !rho baryon.
       thmout(it,5) = bbnvalues%phie          !chemical potential.
       thmout(it,6) = thm(10)       !rho total.
       thmout(it,7) = allnp !n/p ratio
       thmout(it,8)   = bbnvalues%hv/(3.3683e+4_dl)!baryon to photon ratio.
       dtout(it)    = dt            !time step.
       senout(it)   = bbnvalues%sen !entropy ber baryon
       hubout(it)   = hubcst        !expansion rate.

!------equiibrium.

       if (equilflag) call nse_det(tplmev,bbnvalues%y(1),bbnvalues%y(2),thm(9))
       !if (yeflag) call ye_det(tplmev,alln,allp)
       if (yeflag) call ye_det(tcmev,tplmev,alln,allp,yetest3)


#ifdef prllel
       if (rank.eq.0) then
#endif

       sentest = (thm(1) + thm(3) + thm(4) + thm(7))/thm(9)*amu/tplmev
       sentest2 = (thm(1) + thm(3) + thm(4) + thm(7))*(bbnvalues%asf**3)/tplmev
       barnum = thm(9)/mev4tocgs*8._dl*pi/3._dl/mpl**2 &
              *2.99792458_dl**2/197.3269788_dl**2*3.08567758_dl**2*1.e+80_dl*(2.726e-09_dl/bbnvalues%tpl)**3

       write (82,*) tcmev, tplmev, tcmev/tplmev, bbnvalues%asf, dt &
                    ,bbnvalues%sen,thm(8)/thm(1)*8._dl/7._dl*(11._dl/4._dl)**(4._dl/3._dl)
       !write (84,*) t, dt, bbnvalues%tpl, bbnvalues%asf, bbnvalues%phie, bbnvalues%sen, bbnvalues%hv
       write (87,*) tcmev, xout(it,:)
       !write (94,*) tcmev, tplmev, bbnvalues%sen, sentest, sentest2, barnum
       call flush(82)
       !call flush(84)
       call flush(87)
       !write (86,*) tplmev, bbnvalues%y
       !if (nurhoflag) then
         !write (72,*) tcmev, bbnvalues%nuoccprob(1,1)%occfrac(:)
         !write (73,*) tcmev, bbnvalues%nuoccprob(1,2)%occfrac(:)
         if (mswflag) write (74,*) tcmev, bbnvalues%nuoccprob(4,1)%occfrac(:)
         !write (75,*) tcmev, bbnvalues%nuoccprob(4,2)%occfrac(:)
         !call flush(72)
         !call flush(73)
         if (mswflag) call flush(74)
         !call flush(75)
       !end if !nurhoflag

#ifdef prllel
       end if
#endif

       !if (rank.eq.0) write(*,*) tcmev, dt, t


       call bbn_yes_save_state_v2(t,bbnvalues,dt)


       return


       end subroutine bbn_yes_save_out_v2

