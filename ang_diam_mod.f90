
       module ang_diam_mod


!------remarks.

!      Module to control the lowest neutrino mass eigenstate.


!------modules.

       use mainvar
       use sdratiovar


       implicit none


       contains


!-----------------------------------------------------------

       real(dl) function ang_diam(adc,hparm) !Angular Diameter Calculator


!------throughput variables.

       real(dl), intent(in) :: adc !scale factor at photon decoupling
       type(rhovar), intent(in) :: hparm !hubble parameter


!------local variables.

       integer i !index
       integer brnum !Index for Boole's rule
       integer intnum !number of integration points

       real(dl) dasc !width of integration variable
       real(dl) asc !integration variable
       real(dl) intfact !integrand
       real(dl) intsum !integral
       real(dl) hubrate !Hubble rate at a = asc


!------procedure.

       intnum = 1001

       !width of integration variable a:
       dasc = (1._dl - adc)/real(intnum-1, kind=dl)

       intsum = 0._dl
       do i=1,intnum
         brnum = mod(i-2,4) + 1
         !indexed scale factor a:
         asc = real(i-1, kind=dl)/real(intnum-1, kind=dl)*(1._dl - adc) + adc
         !calculate hubrate at a = asc:
         hubrate = hubcalc(asc,hparm)

         intfact = dasc/asc**2/hubrate !units of (Mpc*MeV)^-1

         if (i.eq.1) then
           intfact = brends(1)*intfact
         else if (i.eq.intnum) then
           intfact = brends(2)*intfact
         else
           intfact = brcoeffs(brnum)*intfact
         end if
         intsum = intsum + intfact

       end do !i loop

       ang_diam = intsum*convfact1


       return


       end function ang_diam


       end module ang_diam_mod
