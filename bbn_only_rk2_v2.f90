
       subroutine bbn_only_rk2_v2(dt,y,tpl1,tpl2)
       

!------linkages.

!      called by - [subroutine] bbn_yes_stepper_v2
!      calls     - [subroutine] sol_v2


!------remarks.

!      2nd-order Runge-Kutta method.


       use bbnvar_v2
       use nucsolver


       implicit none


       interface ratenuc_interface
         subroutine ratenuc(t9)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: t9
         end subroutine ratenuc
       end interface ratenuc_interface


       save


!------throughput variables.

       real(dl), intent(in) :: dt !time step
       real(dl), dimension(nnuc), intent(inout) :: y !abundances
       real(dl), intent(in) :: tpl1 !first plasma temp. in 10^9 K
       real(dl), intent(in) :: tpl2 !second plasma temp. in 10^9 K


!------local variables.

       real(dl), dimension(nnuc) :: dydtnuc !abundance derivatives
       real(dl), dimension(nnuc) :: dydt0 !abundance derivatives
       integer i
       real(dl) dnpdt_spectra0
       real(dl) tpl2mev


!------procedure.

       call ratenuc(tpl1)
       call sol_v2(y,dydtnuc,dt,tpl1,1)
       !write (*,*) 'here1',y0(1:2)
       do i=1,size(y)
         y0(i) = y(i)
         dydt0(i) = dydtnuc(i)
         y(i) = y0(i) + dydt0(i)*dt
         if ((y(i).lt.ytmin)) y(i) = ytmin  
       end do


       call ratenuc(tpl2)
       call sol_v2(y,dydtnuc,dt,tpl2,2)
       do i=1,size(y)
         y(i) = y0(i) + 0.5_dl*(dydt0(i) + dydtnuc(i))*dt
         if ((y(i).lt.ytmin)) y(i) = ytmin  
       end do

       dydt(:) = 0.5_dl*(dydt0(:) + dydtnuc(:))
       dydtstep(:) = dydt0(:)


       tpl2mev = boltzmann*tpl2*1.d-06
       !if (heflag.and.(tpl2mev.lt.0.08_dl)) then
       !  y(6) = 1.1_dl*y(6)
       !  heflag = .false.
       !end if

       return


       end subroutine bbn_only_rk2_v2
