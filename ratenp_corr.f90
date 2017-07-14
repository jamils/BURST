
       module ratenp_corr

!------linkages.

!      used by - [module] ratenp
!                [subroutine] bbn_yes_launch_v2


!------remarks.

!      Provides corrections to the n<->p rates.


!------modules.

       use mainvar
       !use gauss
       !use transvar


       implicit none


!------variables.


       contains


!--------------------------------------------------------------------

       real(dl) function get_coul_corr(zarg,xarg,parg)


!------linkages.

!      called by - [functions] urcafc, urcarc, ndecay, indecay
!                              bbn_yes_launch_v2
!      calls     - [functions] gammar, gamcar


!------remarks.

!      Function to compute Coulomb correction to n<->p rates.
!      arXiv: 0905.2781


       implicit none


!------throughput variables.

       real(dl), intent(in) :: zarg !charged lepton sign
       real(dl), intent(in) :: xarg !charged lepton speed
       real(dl), intent(in) :: parg !charged lepton momentum


!------local variables.

       real(dl) sfact
       real(dl) :: cradp = 0.8751_dl/197.3269718_dl
       real(dl) omega, gamf2
       complex(dl) garg1, gamf1


!------procedure.

       sfact = sqrt(1._dl - alphafs**2)
       omega = -zarg/xarg
       garg1 = complex(sfact, omega*alphafs)
       gamf1 = gamcar(garg1)
       gamf2 = gammar(2._dl*sfact + 1._dl)
       get_coul_corr = 2._dl*(1._dl + sfact) &
              *(2._dl*parg*cradp)**(2._dl*(sfact - 1._dl)) &
              !*exp(pi*omega) &
              *exp(pi*omega*alphafs) &
              !*exp(abs(gamf1) - abs(gamf2))
              *exp(2._dl*(abs(gamf1) - abs(gamf2)))

       !get_coul_corr = 2._dl*pi*alphafs/xarg &
       !         /(1._dl - exp(-2._dl*pi*alphafs/xarg))


       !get_coul_corr = 1._dl


       return


       end function get_coul_corr


!--------------------------------------------------------------------

       real(dl) function get_zero_rad_corr(barg,yarg,epsarg)


!------linkages.

!      called by - [functions] urcafc, urcarc, urcafn, urcarn, ndecay, indecay
!                              bbn_yes_launch_v2
!      calls     - none


!------remarks.

!      Function to compute zero temperature radiative correction to n<->p rates.
!      reference: Dicus et al 1982 PRD


       implicit none


!------throughput variables.

       real(dl), intent(in) :: barg !charged lepton speed
       real(dl), intent(in) :: yarg !E_\nu/m_e
       real(dl), intent(in) :: epsarg !E_{e^\pm}/m_e


!------local variables.

       real(dl) rfact
       real(dl) term1, term2, term3 


!------procedure.


       rfact = 1._dl/2._dl/barg*log((1._dl + barg)/(1._dl - barg))
       term1 = 4._dl*(rfact - 1._dl) &
               *(yarg/3._dl/epsarg - 1.5_dl + log(2._dl*yarg))
       term2 = rfact*(2._dl*(1._dl + barg**2) &
               + yarg**2/6._dl/epsarg**2 &
               - 4._dl*barg*rfact)
       term3 = -4._dl*(2._dl + 36._dl*barg + 25._dl*barg**2 &
               + 30._dl*barg**3 + 20._dl*barg**4 + 8._dl*barg**5) &
               /(1._dl + barg)**6

       get_zero_rad_corr = 1._dl + alphafs/2._dl/pi &
                           *(40._dl + term1 + term2 + term3)


       !get_zero_rad_corr = 1._dl


       return


       end function get_zero_rad_corr


       end module ratenp_corr
