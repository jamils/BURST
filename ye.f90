
       module ye


!------modules.

       use mainvar


       implicit none


!------data.


       contains


       subroutine ye_det(tcmev,tpl,yntot,yptot, yetest3)


!------linkages.

!      called by - [subroutine] save_out
!      calls     - none


!-----remarks.

!     Calculates the equilibrium and actual electron fractions.


!------throughput variables.

       real(dl), intent(in) :: tcmev
       real(dl), intent(in) :: tpl !plasma temperature in MeV
       real(dl), intent(in) :: yntot !total neutron abundance
       real(dl), intent(in) :: yptot !total proton abundance
       real(dl), intent(in) :: yetest3


!------data.

       real(dl) :: deltamnp = 1.29333217_dl !mass difference in MeV between proton and neutron


!------local variables.

       real(dl) ye_equil !equilibrium ye
       real(dl) ye_actual !actual ye


!------procedure.

       !equilibrium electron fraction:
       ye_equil = 1._dl/(1._dl + exp(-deltamnp/tpl))

       !actual electron fraction:
       ye_actual = 1._dl/(1._dl + yntot/yptot)

       !write to file:
       !write (23,35) tpl,ye_equil,ye_actual
       write (23,*) tcmev,tpl,ye_equil,ye_actual,yetest3
       call flush(23)

35     format (1p10e12.3)


       return


       end subroutine ye_det


       end module ye
