
       module renorm

!------linkages.

!      used by bbn_yes_therm_v2, ratenp, trans_evolve

!------remarks.

!      Module for calculating renormalization quantities. 

!------modules.

       use mainvar


       contains


!-----------------------------------------------------------

       subroutine get_eg_mass_renorm(tplmev,phie,intscheme,dmeeps2,dmedteps,dmgeps2,dmgdteps)


!------linkages.

!      called by - [subroutine] bbn_no_therm_v2, get_ratenp,
!                               trans_evolve
!      calls     - [none]


!------remarks.

!      computes correction to photon and electron rest masses


!------modules.

!      None


       implicit none


!------throughput variables.

       real(dl), intent(in) :: tplmev !plasma temp.
       real(dl), intent(in) :: phie !electron deg. parameter
       type(bin_scheme), intent(in) :: intscheme
       real(dl), intent(out) :: dmeeps2 !correction to square of xmelec/tplmev^2
       real(dl), intent(out) :: dmedteps !derivative wrt temp of correction to square of xmelec/tplmev^2
       real(dl), intent(out) :: dmgeps2 !correction to square of photon mass/tplmev^2
       real(dl), intent(out) :: dmgdteps !derivative wrt temp of correction to square of photon mass/tplmev^2


!------local variables.

       integer i !indicies
       real(dl) x, dx, intval, intvaldt
       real(dl) mt, occprob1, occprob2


!------procedure division.

!10--------compute factors---------------------------------------


       intval = 0._dl
       intvaldt = 0._dl
       mt = xmelec/tplmev
       do i=1,size(intscheme%abscissas)
         x = intscheme%abscissas(i)
         dx = intscheme%weights(i)
         occprob1 = fd_equil_calc(sqrt(x**2 + mt**2),1._dl,phie)
         occprob2 = fd_equil_calc(sqrt(x**2 + mt**2),1._dl,-phie)

         intval = intval &
                + dx*x**2/sqrt(x**2 + mt**2)*(occprob1 + occprob2)
         intvaldt = intvaldt &
                  + dx*(2._dl*x**2 + mt**2)/sqrt(x**2 + mt**2)*(occprob1 + occprob2)
       end do

       dmeeps2 = 2._dl*pi*alphafs/3._dl &
               + 2._dl*alphafs/pi*intval
       dmedteps = 4._dl*pi*alphafs/3._dl &
                + 2._dl*alphafs/pi*intvaldt

       dmgeps2 = 4._dl*alphafs/pi*intval
       dmgdteps = 4._dl*alphafs/pi*intvaldt


       return


       end subroutine get_eg_mass_renorm


       end module renorm
