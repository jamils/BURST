
       subroutine bbn_yes_driver_v2(t,bbnvalues,dti)


!------linkages.

!      called by - [subroutine] bbn_yes
!      calls     - [subroutine] bbn_yes_stepper_v2, bbn_yes_save_v2


!------remarks.

!      ODE driver for BBN.


       use bbnvar_v2


       implicit none


       interface bbn_yes_stepper_v2_interface
         subroutine bbn_yes_stepper_v2(t,bbnvalues,dt)
           use bbnvar_v2
           implicit none
           real(dl), intent(inout) :: t !time
           type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
           real(dl), intent(inout) :: dt !time step
         end subroutine bbn_yes_stepper_v2
       end interface bbn_yes_stepper_v2_interface

       interface bbn_yes_stepper_old_interface
         subroutine bbn_yes_stepper_old(t,bbnvalues,dt)
           use bbnvar_v2
           implicit none
           real(dl), intent(inout) :: t !time
           type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
           real(dl), intent(inout) :: dt !time step
         end subroutine bbn_yes_stepper_old
       end interface bbn_yes_stepper_old_interface

       interface bbn_yes_save_out_v2_interface
         subroutine bbn_yes_save_out_v2(t,bbnvalues,dt,eqocc)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: t !time
           type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
           real(dl), intent(in) :: dt !time step
           real(dl), dimension(:), intent(in) :: eqocc !equilibrium occupation probabilities
         end subroutine bbn_yes_save_out_v2
       end interface bbn_yes_save_out_v2_interface


       save


!------throughput variables.

       real(dl), intent(inout) :: t !time
       type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
       real(dl), intent(inout) :: dti !initial time step


!------local variables.

       real(dl) tpl, tplmev,tcmev !plasma, comoving temperatures
       real(dl), dimension(:), allocatable :: eqocc !equilibrium occupation probabilities
       real(dl) dt !time step
       real(dl) ndenspro


!------procedure.


!10----input initialization information, relabel


       !dt = 1.e-10_dl
       dt = dti


       tmevi = 1.e-06_dl*boltzmann*bbnvalues%tpl


       allocate(eqocc(size(eps_bins%abscissas))) !allocate eqocc
       eqocc = 1._dl/(exp(eps_bins%abscissas) + 1._dl)


!20----loop until done.

       do !begin looping


!------determine when to turn on transport:

       tplmev = boltzmann*bbnvalues%tpl*1.e-06_dl
       if (tplmev.lt.8._dl) then
         tcmev = (asfdec/bbnvalues%asf)*tcmdec*boltzmann*1.e-06_dl
       else
         tcmev = tplmev
       end if

       nurhoflag = .false.
       if ((transflag).and.(tplmev.lt.8._dl)) nurhoflag = .true.
       if (tcmev.lt.0.015_dl) nurhoflag = .false.


!------accumulate.

       tpl = bbnvalues%tpl
       !if ((tpl.le.t9f).or.(dt.lt.abs(dtlow/dlt9dt))) then
       if (tpl.le.t9f) then
         call bbn_yes_save_out_v2(t,bbnvalues,dt,eqocc)
         exit
       end if

       if (ip.eq.inc) then
         call bbn_yes_save_out_v2(t,bbnvalues,dt,eqocc)
         if (it.eq.itmax) exit
       end if                 


!------compute derivatives of variables to be evolved.

       call bbn_yes_stepper_v2(t,bbnvalues,dt)


!------reset counters.

       ndenspro = sum(bbnvalues%y(:)*zm(:))*thm(9)/amu/mev4tocgs
#ifdef prllel
       if (rank.eq.0) then
#endif
       !write (95,*) tcmev, tplmev, ndensem, ndensep, ndenspro, bbnvalues%phie &
       !   ,(ndensem - ndensep - ndenspro)/(ndensem + ndensep + ndenspro)
       !write (96,*) tcmev, tplmev, bbnvalues%sen, dspldt &
       !               ,snu, dsnudt, stot, dstotdt
       !call flush(95)
       !call flush(96)
       !if (it.ge.249) write (*,*) 'here1', dt, ip, it
#ifdef prllel
       end if
#endif

       if (ip.eq.inc) ip = 0 !reset iteration counters.
       ip = ip + 1
       is = is + 1

       !if (it.eq.50) stop

       end do !do loop without index


       bbnvalues%nuoccprob(:,:)%dilfact = tcmev/tplmev

       deallocate(eqocc)

       
       return


       end subroutine bbn_yes_driver_v2

