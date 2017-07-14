
       subroutine bbn_yes_rkck_v2(t,bbntemp,dt,bbnerr,bbnscal)
       

!------linkages.

!      called by - [subroutine] bbn_yes_stepper_v2
!      calls     - [subroutine] bbn_yes_evolve_v2


!------remarks.

!      Explicit fifth-order Runge-Kutta algorithm.
!      Uses Cash-Karp method to take a step, and
!      embedded fourth-order method for error monitoring.


       use bbnvar_v2


       implicit none


       save


       interface bbn_yes_evolve_v2_interface
         subroutine bbn_yes_evolve_v2(t,bbnvalues,bbnderivs,loop)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: t !time
           type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution values
           type(bbnevolvar), intent(out) :: bbnderivs !bbn evolution variable derivatives
           integer, intent(in) :: loop
         end subroutine bbn_yes_evolve_v2
       end interface bbn_yes_evolve_v2_interface



!------throughput variables.

       real(dl), intent(inout) :: t !time
       type(bbnevolvar), intent(inout) :: bbntemp !bbn evolution values
       real(dl), intent(in) :: dt !time step
       type(bbnevolvar), intent(out) :: bbnerr !error in bbn evolution values
       type(bbnevolvar), intent(out) :: bbnscal !error in bbn evolution values


!------local variables.

       real(dl) tin
       !type(bbnevolvar), pointer :: bbnvaluesin
       type(bbnevolvar) :: bbnvaluesin
       !type(bbnevolvar), dimension(:), pointer :: bbnderivs !array of derivatives
       type(bbnevolvar), dimension(6) :: bbnderivs !array of derivatives

       integer i, j, l, m, n !indices

       real(dl), dimension(6) :: addx = (/0._dl, 0.2_dl, 0.3_dl, 0.6_dl, 1._dl, 0.875_dl/)
       real(dl), dimension(5,5) :: inty = transpose(reshape([ &
               0.2_dl,             0._dl,           0._dl,             0._dl,                0._dl, &
               3._dl/40._dl,       9._dl/40._dl,    0._dl,             0._dl,                0._dl, &
               0.3_dl,             -0.9_dl,         1.2_dl,            0._dl,                0._dl, &
               -11._dl/54._dl,     2.5_dl,          -70._dl/27._dl,    35._dl/27._dl,        0._dl, &
               1631._dl/55296._dl, 175._dl/512._dl, 575._dl/13824._dl, 44275._dl/110592._dl, 253._dl/4096._dl], &
               shape(inty)))
       real(dl), dimension(6), parameter :: accumy = (/37._dl/378._dl, 0._dl, 250._dl/621._dl, &
                                                       125._dl/594._dl, 0._dl, 512._dl/1771._dl/)
       real(dl), dimension(6) :: diff54 = (/accumy(1) - 2825._dl/27648._dl, 0._dl, &
                                            accumy(3) - 18575._dl/48384._dl, accumy(4) - 13525._dl/55296._dl, &
                                            -277._dl/14336._dl, accumy(6) - 0.25_dl/)


!------procedure.

       !if (.not.allocated(bbnvaluesin%nuoccprob)) then
       if (.not.associated(bbnvaluesin%nuoccprob)) then
         allocate(bbnvaluesin%nuoccprob(size(bbntemp%nuoccprob,1),2))
         do i=1,6
           allocate(bbnderivs(i)%nuoccprob(size(bbntemp%nuoccprob,1),2))
         end do !i
         do m=1,size(bbntemp%nuoccprob,1)
           do n=1,2
             allocate(bbnvaluesin%nuoccprob(m,n)%occfrac(size(bbntemp%nuoccprob(1,1)%occfrac)))
             do i=1,6
               allocate(bbnderivs(i)%nuoccprob(m,n)%occfrac(size(bbntemp%nuoccprob(1,1)%occfrac)))
             end do !i
           end do !n
         end do !m
       end if


       do l=1,6
         tin = t + addx(l)*dt
         call bbn_yes_assigntype_v2(bbnvaluesin,bbntemp)
         if (l.ge.2) then
           do m=1,(l-1)
             call bbn_yes_addtype_v2(bbnvaluesin,bbnderivs(m),inty(l-1,m)*dt)
             !if (isnan(bbnvaluesin%tpl)) write(*,*) is, l
           end do
         end if
         call bbn_yes_evolve_v2(tin,bbnvaluesin,bbnderivs(l),l)
       end do


       do l=1,6
         call bbn_yes_addtype_v2(bbntemp,bbnderivs(l),accumy(l)*dt)
         call bbn_yes_addtype_v2(bbnerr,bbnderivs(l),diff54(l)*dt)
       end do
       
       !do m=1,size(bbnscal%nuoccprob, 1)
       !  do n=1,2
       !    do i=1,size(bbnscal%nuoccprob(1,1)%occfrac)
       !      bbnscal%nuoccprob(m,n)%occfrac(i) = abs(dt*bbnderivs(1)%nuoccprob(m,n)%occfrac(i))
             !bbnscal%nuoccprob(m,n)%occfrac(i) = 1.e+10_dl
       !    end do !i
       !  end do !n
       !end do !m
       
       
       return


       end subroutine bbn_yes_rkck_v2

