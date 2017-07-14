
       module opt_depth_mod


!------remarks.

!      Module to calculate the optical depth for given xec.


!------modules.

       use mainvar
       use sdratiovar


       implicit none


       contains


!--------------------------------------------------------------------

       subroutine find_adc(ascarray,xec,hparm,adcod,atest)


!------linkages.

!      called by - [subroutine] recom


!------remarks.

!      Computes recombination history over a specified epoch.


       implicit none


       save


!------throughput variables.

       real(dl), dimension(:), intent(in) :: ascarray !Scale factor values for xec
       real(dl), dimension(:), intent(in) :: xec !Free-electron fraction
       type(rhovar), intent(in) :: hparm !hubble parameter
       real(dl), intent(out) :: adcod !scale factor at decoupling
       real(dl), intent(inout), optional :: atest !guess for adc


!------local variables.

       real(dl) aguess !guess for adc
       real(dl) alow !lower estimate for adc
       real(dl) ahigh !upper estimate for adc
       real(dl) ain !scale factor to iterate on
       real(dl) tau_od !optical depth
       integer counter !used in case integrals do not converge

!------procedure.

       if (.not.present(atest)) then
         aguess = 9.15e-04_dl
       else
         aguess = atest
       end if

       alow = 0.9_dl*aguess
       counter = 1
       do
         tau_od = opt_depth(alow,ascarray,xec,hparm)
         if ((tau_od.gt.1._dl).or.(counter.eq.100)) then
           exit
         else
           alow = 0.9_dl*alow
         end if
         counter = counter + 1
       end do

       ahigh = 1.1_dl*aguess
       counter = 1
       do
         tau_od = opt_depth(ahigh,ascarray,xec,hparm)
         if ((tau_od.lt.1._dl).or.(counter.eq.100)) then
           exit
         else
           ahigh = 1.1_dl*ahigh
         end if
         counter = counter + 1
       end do

       ain = aguess
       counter = 1
       do
         tau_od = opt_depth(ain,ascarray,xec,hparm)
         if ((abs(tau_od-1._dl).lt.1.e-06_dl).or.(counter.eq.100)) then
           exit
         else
           if (tau_od.gt.1._dl) then
             alow = ain
             ain = (ain + ahigh)/2._dl
           else
             ahigh = ain
             ain = (ain + alow)/2._dl
           end if
         end if
         counter = counter + 1
       end do

       adcod = ain

#ifdef prllel
       if (rank.eq.0) then
#endif

       !write (*,*) 'a_{dc} = ', ain
       !write (*,*) 'Optical Depth:', tau_od
       !write (*,*) 'counter:', counter

#ifdef prllel
       end if !(rank.eq.0)
#endif


       return


       end subroutine find_adc


!-----------------------------------------------------------

       real(dl) function opt_depth(ain,ascarray,xec,hparm) !Optical Depth Calculator


!------modules.

       use interp
       use xe_history


!------throughput variables.

       real(dl), intent(in) :: ain !input scale factor
       real(dl), dimension(:), intent(in) :: ascarray !Scale factor values for xec
       real(dl), dimension(:), intent(in) :: xec !Free-electron fraction
       type(rhovar), intent(in) :: hparm !hubble parameter


!------local variables.

       integer i !index
       integer brnum !Index for Boole's rule
       integer intnumod !number of integration points

       real(dl) dasc !width of integration variable
       real(dl) asc !integration variable
       real(dl) intfact !integrand
       real(dl) intsum !integral
       real(dl) hubrate !Hubble rate at a = asc

       real(dl), dimension(:), allocatable :: xecddot !Second derivatives of xec
       real(dl), dimension(3) :: xeitemp !temp. array of xei at single value of scale factor
       real(dl) xetot !Total ionization fraction


!------procedure.

       !allocate(xecddot(size(xec))) !Declare size of xecddot
       allocate(xecddot(lastind)) !Declare size of xecddot

       intnumod = 10001

       !width of integration variable a:
       dasc = (1._dl - ain)/real(intnumod-1, kind=dl)

       call cubic_spline(ascarray(1:lastind),xec(1:lastind),xecddot,0._dl,0._dl)
       !do i=1,lastind
       !  write (101,*) ascarray(i), xec(i), xecddot(i)
       !end do
       !call flush(101)

       intsum = 0._dl
       do i=1,intnumod

         brnum = mod(i-2,4) + 1
         !indexed scale factor a:
         asc = real(i-1, kind=dl)/real(intnumod-1, kind=dl)*(1._dl - ain) + ain
         !calculate hubrate at a = asc:
         hubrate = hubcalc(asc,hparm)

         !if (asc.gt.ascarray(size(ascarray))) then
         if (asc.gt.ascarray(lastind)) then
           !xetot = xec(size(xec))
           xetot = xec(lastind)
         else if (asc.lt.ascarray(1)) then
           call saha_net(asc,xeitemp)
           xetot = sum(xeitemp*xeiz)
         else
           !call cubic_splint(ascarray,xec,xecddot,asc,xetot)
           !call cubic_splint(ascarray(1:lastind),xec(1:lastind),xecddot,asc,xetot)
           xetot = cubic_splint(ascarray(1:lastind),xec(1:lastind),xecddot,asc)
         end if
         !do i=1,lastind
         !  write (101,*) asc, xetot
         !end do
         call flush(101)

         intfact = ne0*(1._dl/asc)**3*xetot*sigmat*dasc/asc/hubrate !dimensionless

         if (i.eq.1) then
           intfact = brends(1)*intfact
         else if (i.eq.intnumod) then
           intfact = brends(2)*intfact
         else
           intfact = brcoeffs(brnum)*intfact
         end if
         intsum = intsum + intfact

       end do !i loop


       opt_depth = intsum


       deallocate(xecddot)


       return


       end function opt_depth


       end module opt_depth_mod
