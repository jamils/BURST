
       module nu_assign


!------linkages.

!      used by - [program] main


!------remarks.

!      Assigns occupation probabilities


!------modules.

       use mainvar


       implicit none


       contains


!--------------------------------------------------------------------

       subroutine nu_assign_equil(nuarray,xiin,mass,tempratioin)


!------linkages.

!      called by - [program] main
!      calls     - none

!------remarks.

!      Subroutine to assign values to a type(nuvar) variable.

       use mainvar


       implicit none


!------throughput variables.

       type(nuvar), intent(out) :: nuarray !nu parameters
       real(dl), intent(in) :: xiin !input neutrino degeneracy parameter
       real(dl), intent(in) :: mass !mass in eV
       real(dl), intent(in) :: tempratioin !Parameter used for temperature ratio


!------local variables.

       integer i !indicies
       real(dl) tempratio


!------procedure.


       !nbins = size(eps_bins)


!------nu occupation fractions.

       !nuarray%occfrac = 1._dl/(exp(eps_bins%abscissas - xiin) + 1._dl) !spectrum

       nuarray%mass = mass !mass
       nuarray%dilfact = 1._dl !dilution factor
       nuarray%mix = 1._dl !\sin^2\theta

       if ((tempratioin.le.0._dl).or.(tempratioin.gt.1._dl)) then
         !18Apr2015: put in write statement to error log.
         tempratio = 1._dl
       else
         tempratio = tempratioin
       end if

       do i=1,nbins
         nuarray%occfrac(i) = 1._dl/(exp(eps_bins%abscissas(i)/tempratio - xiin) + 1._dl) !spectrum
       end do


       return


       end subroutine nu_assign_equil


!--------------------------------------------------------------------

       subroutine nu_assign_null(nuarray,mass,s2mixangle)


!------linkages.

!      called by - [program] main
!      calls     - none

!------remarks.

!      Subroutine to assign null values to a type(nuvar) variable.

       use mainvar


       implicit none


!------throughput variables.

       type(nuvar), intent(out) :: nuarray !nu parameters
       real(dl), intent(in) :: mass !mass in eV
       real(dl), intent(in) :: s2mixangle !\sin(2\theta)


!------local variables.


!------procedure.


       !nbins = size(eps_bins)


!------nu occupation fractions.

       !nuarray%occfrac = 1._dl/(exp(eps_bins%abscissas - xiin) + 1._dl) !spectrum

       nuarray%mass = mass !mass
       nuarray%dilfact = 1._dl !dilution factor
       nuarray%mix = s2mixangle !\sin(2\theta)

       nuarray%occfrac = 0._dl


       return


       end subroutine nu_assign_null


       end module nu_assign
