
       subroutine main_init


!------linkages.

!      called by - [program] main
!      calls     - [subroutine] bbn_init, trans_init, sdratio_init


!------remarks.

!      Initiates the variables for BBN, trans, ratio, ....


!------modules.

       !use bbnvar
       !use sdratiovar


       implicit none

       !interface trans_interface
       !  subroutine trans_init(nutrans)
       !    use transvar
       !    implicit none
       !    type(nuvar), dimension(:,:), intent(inout) :: nutrans !neutrino information for trans
       !  end subroutine trans_init
       !end interface trans_interface

!------local variables.


!------procedure.

       call bbn_init
       call trans_init
       call sdratio_init

       !Put future bsm physics (...)_init calls here, e.g.:
       !call darkmatter_init


       return


       end subroutine main_init
