
       subroutine sdratio_launch(obh2in,neffparm,ypri,tneff)
       
!------linkages.

!      called by - [subroutine] sdratio
!      calls     - [subroutine] integ_a
!                  [function] rhonulagint


!-----remarks.

!     Computes the sound horizon and diffusion length at recombination,
!     relating the two quantities to Neff


!------modules.

       use sdratiovar
       use sdratio_length
       !use interp
       !use massvar
       use xe_history
       !use ang_diam_mod


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: obh2in !input baryon number
       real(dl), intent(in) :: neffparm !parameter for Neff
       real(dl), intent(in) :: ypri !primordial 4He mass fraction
       real(dl), intent(inout) :: tneff !Neff target in; Neff calculated out


!------local variables.

       real(dl) drhorad0




!------procedure.

       obh2 = obh2in !for book keeping purposes
       !omegag = hubparm%og

       neffrangel = tneff - 0.5_dl*neffwidth
       if (neffrangel.lt.0._dl) neffrangel = 0._dl
       neffrangeh = tneff + 0.5_dl*neffwidth

       !number density of baryons at current epoch
       !in units of MeV**3:
       nb0 = obh2/h**2*hubparm%rhoc0/amu
       !energy density of baryons at current epoch
       !in units of MeV**4:
       rhob0 = amu*nb0

       !och2 = och2 + 0.01_dl !only used for testing purposes
       hubparm%om = (och2 + obh2)/h**2
       !omegam = (och2 + obh2)/h**2

       ypxe = ypri !Y_p needed for recombination history.
       !number density of all electrons at current epoch:
       ne0 = nb0*(1._dl - ypxe/2._dl)
       !Max values for ratio of ions to total electrons
       xeimax(1) = (1._dl - ypxe)/(1._dl - ypxe/2._dl) !HII
       xeimax(2) = ypxe/4._dl/(1._dl - ypxe/2._dl) !HeII
       xeimax(3) = ypxe/4._dl/(1._dl - ypxe/2._dl) !HeIII
       xeimax(4) = 0._dl !baryon temp.

       !atranshe32 = sinrsolver(saha_ion,dsida,1.e-04_dl*xeimax(3) &
       !        ,1.e-06_dl,2.e-04_dl,ypxe/2._dl,1._dl,ipothe2,xeimax(3))
       atranshe32 = sinrsolver(saha_ion,dsida,1.e-04_dl*xeimax(3) &
               ,1.e-06_dl,2.e-04_dl,1._dl,ipothe2,xeimax(3),xeimax(3))
       !atranshe21 = sinrsolver(saha_ion,dsida,0.99_dl*xeimax(2) &
       !        ,1.e-06_dl,4.e-04_dl,ypxe/2._dl,4._dl,ipothe1,xeimax(2))
       atranshe21 = sinrsolver(saha_ion,dsida,0.99_dl*xeimax(2) &
               ,1.e-06_dl,4.e-04_dl,4._dl,ipothe1,2._dl*xeimax(2),xeimax(2))
       !most likely won't be able to use this:
       !atransh21 = sinrsolver(saha_ion,dsida,0.99_dl*xeimax(1) &
       !        ,1.e-06_dl,6.e-04_dl,(1._dl - ypxe/2._dl),1._dl,ipoth1,xeimax(1))
       atransh21 = sinrsolver(saha_ion,dsida,0.99_dl*xeimax(1) &
               ,1.e-06_dl,6.e-04_dl,1._dl,ipoth1,(xeimax(1) + 2._dl*xeimax(2) - 0.016_dl),xeimax(1))


       !hubparm%orad = 0._dl
       if (neffparm.ge.0._dl) then
         drhorad0 = neffparm*t0mev**4
       else
         drhorad0 = 0._dl
       end if
       !omegarad = drhorad0/hubparm%rhoc0
       hubparm%orad = drhorad0/hubparm%rhoc0

       xei => ratioy(:,1:3)
       tbar => ratioy(:,4)

 
       return


       end subroutine sdratio_launch
