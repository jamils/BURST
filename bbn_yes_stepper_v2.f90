
       subroutine bbn_yes_stepper_v2(t,bbnvalues,dt)


!------linkages.

!      called by - [subroutine] bbn_yes_driver_v2
!      calls     - [subroutine] bbn_yes_rkck_v2


!------remarks.

!      Executes a single step through BBN evolution.


       use bbnvar_v2
       use msw


       implicit none


       save


       interface bbn_yes_rkck_v2_interface
         subroutine bbn_yes_rkck_v2(t,bbntemp,dt,bbnerr,bbnscal)
           use bbnvar_v2
           implicit none
           real(dl), intent(inout) :: t !time
           type(bbnevolvar), intent(inout) :: bbntemp !bbn evolution values
           real(dl), intent(in) :: dt !time step
           type(bbnevolvar), intent(out) :: bbnerr !error in bbn evolution values
           type(bbnevolvar), intent(out) :: bbnscal !error in bbn evolution values
         end subroutine bbn_yes_rkck_v2
       end interface bbn_yes_rkck_v2_interface

       interface bbn_only_rk2_v2_interface
         subroutine bbn_only_rk2_v2(dt,y,tpl1,tpl2)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: dt !time step
           real(dl), dimension(nnuc), intent(inout) :: y !abundances
           real(dl), intent(in) :: tpl1 !first plasma temp. in 10^9 K
           real(dl), intent(in) :: tpl2 !second plasma temp. in 10^9 K
         end subroutine bbn_only_rk2_v2
       end interface bbn_only_rk2_v2_interface

!------throughput variables.

       real(dl), intent(inout) :: t !time
       type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
       real(dl), intent(inout) :: dt !time step


!------local variables.

       type(bbnevolvar) :: bbnscal
       type(bbnevolvar) :: bbntemp
       type(bbnevolvar) :: bbnerr
       real(dl) dt0, dttemp !time step variables
       real(dl) tol, errmax !quality-control variables
       real(dl) err, scal !shorthand for bbnerr%nuoccprob%occfrac and bbnscal
       integer i,m,n !indices
       real(dl) dtmin, dtl
       real(dl) dtold
       real(dl) tplmev
       integer dttrip,dttrip2
       real(dl) dtmsw !time step from msw


!------procedure.

       if (.not.associated(bbnscal%nuoccprob)) then
         allocate(bbnscal%nuoccprob(size(bbnvalues%nuoccprob,1),2)) !allocate bbnscal
         allocate(bbntemp%nuoccprob(size(bbnvalues%nuoccprob,1),2)) !allocate bbntemp
         allocate(bbnerr%nuoccprob(size(bbnvalues%nuoccprob,1),2)) !allocate bbnerr
         do m=1,size(bbnvalues%nuoccprob,1)
           do n=1,2
             allocate(bbnscal%nuoccprob(m,n)%occfrac(size(bbnvalues%nuoccprob(1,1)%occfrac)))
             allocate(bbntemp%nuoccprob(m,n)%occfrac(size(bbnvalues%nuoccprob(1,1)%occfrac)))
             allocate(bbnerr%nuoccprob(m,n)%occfrac(size(bbnvalues%nuoccprob(1,1)%occfrac)))
           end do !n
         end do !m
       end if

       tol = 1.e-08_dl

       dt0 = dt

       call bbn_yes_assigntype_v2(bbnscal,bbnvalues)

       bbnscal%phie = 17._dl
       bbnerr%y = 0._dl


       call bbn_yes_assigntype_v2(bbntemp,bbnvalues)
       do
         bbnerr%tpl = 0._dl
         bbnerr%hv = 0._dl
         bbnerr%phie = 0._dl
         bbnerr%asf = 0._dl
         bbnerr%sen = 0._dl
         do m=1,size(bbnvalues%nuoccprob,1)
           do n=1,2
             do i=1,size(bbnvalues%nuoccprob(1,1)%occfrac)
               bbnerr%nuoccprob(m,n)%occfrac(i) = 0._dl
             end do !i
           end do !n
         end do !m
         call bbn_yes_rkck_v2(t,bbntemp,dt,bbnerr,bbnscal)
         errmax = 0._dl
         errmax = max(errmax,abs(bbnerr%tpl/bbnscal%tpl))
         errmax = max(errmax,abs(bbnerr%hv/bbnscal%hv))
         errmax = max(errmax,abs(bbnerr%phie/bbnscal%phie))
         errmax = max(errmax,abs(bbnerr%asf/bbnscal%asf))
         errmax = max(errmax,abs(bbnerr%sen/bbnscal%sen))
         dttrip = 0
         tplmev = bbntemp%tpl*boltzmann*1.e-06_dl
         !if (tplmev.gt.0.05_dl) then
           do m=1,size(bbnvalues%nuoccprob,1)
             do n=1,2
               do i=1,size(bbnvalues%nuoccprob(1,1)%occfrac)
                 err = bbnerr%nuoccprob(m,n)%occfrac(i)
                 scal = bbnscal%nuoccprob(m,n)%occfrac(i)
                 errmax = max(errmax,abs(err/scal))
               end do !i
             end do !n
           end do !m
         !end if !tplmev
         errmax = errmax/tol
         if (errmax.gt.1._dl) then !did not meet error tolerance; redo
           call bbn_yes_assigntype_v2(bbntemp,bbnvalues)
           dttemp = 0.9_dl*dt/errmax**0.25_dl
           !dtmin = abs(1._dl/dlt9dt)*ct  
           !dt = min(dt, dtmin)
           dt = max(dttemp, 0.1_dl*dt)
           if (dt.lt.(1.e-08_dl*dt0)) then
             write (*,*) 'issue in bbn_stepper'
             exit
           end if
         else !met error tolerance
           yetest3 = yetest3 +dyetest3dt*dt
           snu = snu + dsnudt*dt
           stot = stot + dstotdt*dt
           !insert call to nucsolver
           call bbn_only_rk2_v2(dt,bbntemp%y,bbnvalues%tpl,bbntemp%tpl)

           if (mswflag) then
             call msw_main(bbntemp%nuoccprob &
                  ,tcmmsw,bbnvalues%tpl*boltzmann*1.e-06_dl &
                  ,yemsw,rhoemsw,nepdiffmsw,nbarymsw &
                  ,dtpldtmsw,dyedtmsw,drhoedtmsw,dnepdtmsw &
                  ,hubmsw,dtmsw)
           end if

           t = t + dt
           call bbn_yes_assigntype_v2(bbnvalues,bbntemp)
           dtold = dt
           if (errmax.gt.1.89e-04_dl) then
             dt = 0.9_dl*dt/errmax**0.2_dl
           else
             if (is.gt.3) then
               dt = 5._dl*dt
             else
               dt = abs(1._dl/dlt9dt)*ct
             end if
           end if
           dtmin = dt
           tplmev = bbnvalues%tpl*boltzmann*1.e-06_dl
           if (is.gt.3) then
           !if (tplmev.lt.10._dl) then
             do i=1,isize
               if ((dydtstep(i).ne.0._dl).and.(bbntemp%y(i).gt.ytmin)) then
                 dtl = abs(bbntemp%y(i)/dydtstep(i))*cy &
                       *(1._dl+(dlog10(bbntemp%y(i))/dlog10(ytmin))**2)  !(ref 2).
                 !dtl = abs(bbntemp%y(i)/dydtstep(i))*cy
                 !if (dtl.lt.dtmin) dtmin = dtl
                 if (dtl.lt.dtmin) then
                   dtmin = dtl
                   dttrip2 = i
                 end if
               end if
             end do
           end if
           !if (.not.mswflag) then
             dt = min(dtmin,dt)
           !else
           !  dt = min(dtmin,dt,dtmsw)
           !end if
           !if (dt.lt.1.e-03_dl*t) then
           !  dt = 1.e-03_dl*t
           !end if
           exit
         end if
       end do



       !deallocate(bbnscal%nuoccprob)
       !deallocate(bbntemp%nuoccprob)
       !deallocate(bbnerr%nuoccprob)
       !nullify(bbnscal%nuoccprob)
       !nullify(bbntemp%nuoccprob)
       !nullify(bbnerr%nuoccprob)


       return


       end subroutine bbn_yes_stepper_v2

