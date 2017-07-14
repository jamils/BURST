
       subroutine bbn_yes_therm_v2(bbnvalues, tcmev, loop)

!------linkages.

!      called by - [subroutine] bbn_yes_evolve_v2
!      calls     - [subroutine] bessel_eval


!------remarks.

!      computes various temperature dependent thermodynamic quantities,
!      without BBN.


!------modules.

       use bbnvar_v2
       use bessel
       use renorm


       implicit none


       !interface bbn_no_renorm_v2_interface
       !  subroutine bbn_no_renorm_v2(tplmev,phie,intscheme,dmeeps2,dmedteps,dmgeps2,dmgdteps)
       !    use mainvar
       !    implicit none
       !    real(dl), intent(in) :: tplmev !plasma temp.
       !    real(dl), intent(in) :: phie !electron deg. parameter
       !    type(bin_scheme), intent(in) :: intscheme
       !    real(dl), intent(out) :: dmeeps2 !correction to square of xmelec/tplmev^2
       !    real(dl), intent(out) :: dmedteps !derivative wrt temp of correction to square of xmelec/tplmev^2
       !    real(dl), intent(out) :: dmgeps2 !correction to square of photon mass/tplmev^2
       !    real(dl), intent(out) :: dmgdteps !derivative wrt temp of correction to square of photon mass/tplmev^2
       !  end subroutine bbn_no_renorm_v2
       !end interface bbn_no_renorm_v2_interface


!------throughput variables.

       type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
       real(dl), intent(in) :: tcmev !comoving temperature in MeV
       integer, intent(in) :: loop !loop number


!------local variables.

       integer i,m,n !indicies
       real(dl) tpl, tplmev !plasma temperatures
       real(dl) asf, phie, hv
       real(dl) z !defined by z = m(electron)*c**2/k*tpl.
       real(dl) rhonu !\nu relativistic energy density
       real(dl) rhonui !\nu_i relativistic energy density
       real(dl) rhor !Extra relativistic energy density
       real(dl) cosh1,cosh2,cosh3,cosh4,cosh5
       real(dl) sinh1,sinh2,sinh3,sinh4,sinh5

       real(dl) glag1, glag2, glag3
       real(dl) glag4, glag5, glag6, glag7
       real(dl) glag11, glag12
       real(dl) mt, x, dx, occprob
       real(dl) occprob1, occprob2
       real(dl) dmeeps2, dmedteps
       real(dl) dmgeps2, dmgdteps
       real(dl) glagnem, glagnep


!------procedure division.

!10--------compute factors---------------------------------------

       tpl = bbnvalues%tpl
       hv = bbnvalues%hv
       !z = 5.930_dl/tpl !z = m(electron)c**2/k(tpl).
       tplmev = boltzmann*tpl*1.e-06_dl
       phie = bbnvalues%phie
       asf = bbnvalues%asf



!..........trignometric function values.
       if (phie.le.17._dl) then        !no chance of overflow.
         cosh1 = cosh(phie)
         cosh2 = cosh(2._dl*phie)
         cosh3 = cosh(3._dl*phie)
         cosh4 = cosh(4._dl*phie)
         cosh5 = cosh(5._dl*phie)
         sinh1 = sinh(phie)
         sinh2 = sinh(2._dl*phie)
         sinh3 = sinh(3._dl*phie)
         sinh4 = sinh(4._dl*phie)
         sinh5 = sinh(5._dl*phie)
       else
         cosh1 = 0._dl
         cosh2 = 0._dl
         cosh3 = 0._dl
         cosh4 = 0._dl
         cosh5 = 0._dl
         sinh1 = 0._dl
         sinh2 = 0._dl
         sinh3 = 0._dl
         sinh4 = 0._dl
         sinh5 = 0._dl
       end if       

!20--------compute thermodynamic variables--------------------

       !thm(1)  = 8.418_dl*tpl**4 !(ref 1).
       !thm(1)  = mev4tocgs*pi**2/15._dl*tplmev**4 !(ref 1).
       !thm(2)  = 4._dl*thm(1)/tpl !(ref 2).
       !thm(3)  = thm(1)/3._dl !(ref 3).

       !thm(4)  = 3206._dl*(bmz(1) - bmz(2) + bmz(3) &
       !          - bmz(4) + bmz(5))                !(ref 4).
       !thm(5)  = 3206._dl*(z/tpl)*(bnz(1) - 2._dl*bnz(2) &
       !          + 3._dl*bnz(3) - 4._dl*bnz(4) + 5._dl*bnz(5)) !(ref 5).
       !thm(6)  = 3206._dl*(bmz(1)*sinh1 - 2._dl*bmz(2)*sinh2 + 3._dl*bmz(3)*sinh3 &
       !          - 4._dl*bmz(4)*sinh4 + 5._dl*bmz(5)*sinh5)                !(ref 6).
       !thm(7)  = 3206._dl*(blz(1)/z - blz(2)/(2._dl*z) &
       !          + blz(3)/(3._dl*z) - blz(4)/(4._dl*z) &
       !          + blz(5)/(5._dl*z))                      !(ref 7).

       rhonu = 0._dl
       do m=1,size(bbnvalues%nuoccprob, 1)
         do n=1,2
           rhonu = rhonu + 1._dl/2._dl/pi**2*tcmev**4 & !units of MeV^4
                 *sum(eps_bins%weights(:)*eps_bins%abscissas(:)**3 &
                 *bbnvalues%nuoccprob(m,n)%occfrac(:))
         end do
       end do
       !rhonu = 3._dl*7._dl/8._dl*pi**2/15._dl*tcmev**4
       thm(8) = mev4tocgs*rhonu !units of g/cm^3

       !thm(9)  = hv*tpl**3
       thm(9)  = rhob0/asf**3
       !(ref 9).

       !rhor = rhors*rnb**(4._dl/3._dl)
       rhor = rhors/asf**4


       thm(13) = 0._dl
       thm(14) = 0._dl

       !call get_eg_mass_renorm(tplmev,phie,glagbessel,dmeeps2,dmedteps,dmgeps2,dmgdteps)
       dmeeps2 = 0._dl
       dmedteps = 0._dl
       dmgeps2 = 0._dl
       dmgdteps = 0._dl

       glag1 = 0._dl
       glag2 = 0._dl
       glag3 = 0._dl
       glag4 = 0._dl
       glag5 = 0._dl
       glag6 = 0._dl
       glag7 = 0._dl
       glag11 = 0._dl
       glag12 = 0._dl
       glagnem = 0._dl
       glagnep = 0._dl
       mt = xmelec/tplmev
       do i=1,size(glagbessel%abscissas)
         x = glagbessel%abscissas(i)
         dx = glagbessel%weights(i)

         occprob = be_equil_calc(sqrt(x**2 + dmgeps2),1._dl,0._dl)

         glag1 = glag1 &
               + dx*x**2*sqrt(x**2 + dmgeps2)*occprob
         glag2 = glag2 &
               + dx*sqrt(x**2 + dmgeps2)*(4._dl*x**2 + dmgeps2)*occprob &
               - dx*dmgdteps*sqrt(x**2 + dmgeps2)*occprob
         glag3 = glag3 &
               + dx*x**4/3._dl/sqrt(x**2 + dmgeps2)*occprob

         occprob1 = fd_equil_calc(sqrt(x**2 + mt**2 + dmeeps2),1._dl,phie)
         occprob2 = fd_equil_calc(sqrt(x**2 + mt**2 + dmeeps2),1._dl,-phie)

         glag4 = glag4 &
               + dx*x**2*sqrt(x**2 + mt**2 + dmeeps2)*(occprob1 + occprob2)
         glag5 = glag5 &
               + dx*sqrt(x**2 + mt**2 + dmeeps2)*(4._dl*x**2 + mt**2 + dmeeps2)*(occprob1 + occprob2) &
               - dx*dmedteps*sqrt(x**2 + mt**2 + dmeeps2)*(occprob1 + occprob2)
         glag6 = glag6 &
               + dx*(3._dl*x**2 + mt**2 + dmeeps2)*(occprob1 - occprob2)
         glag7 = glag7 &
               + dx*x**4/3._dl/sqrt(x**2 + mt**2 + dmeeps2)*(occprob1 + occprob2)
         glag11 = glag11 &
               + dx*(mt**2 + dmeeps2)*(occprob1 - occprob2) &
               - dx*dmedteps/2._dl*(occprob1 - occprob2)
         glag12 = glag12 &
               + dx*(2._dl*x**2 + mt**2 + dmeeps2)/sqrt(x**2 + mt**2 + dmeeps2)*(occprob1 + occprob2)

         glagnem = glagnem &
                 + dx*x**2*occprob1
         glagnep = glagnep &
                 + dx*x**2*occprob2
       end do
       glag1 = 1._dl/pi**2*tplmev**4*glag1*mev4tocgs
       glag2 = 1._dl/pi**2*tplmev**3*glag2*mev3tocgs
       glag3 = 1._dl/pi**2*tplmev**4*glag3*mev4tocgs
       glag4 = 1._dl/pi**2*tplmev**4*glag4*mev4tocgs
       glag5 = 1._dl/pi**2*tplmev**3*glag5*mev3tocgs
       glag6 = 1._dl/pi**2*tplmev**4*glag6*mev4tocgs
       glag7 = 1._dl/pi**2*tplmev**4*glag7*mev4tocgs
       glag11 = 1._dl/tpl/2._dl*glag11
       glag12 = 1._dl/2._dl*glag12
       ndensem = 1._dl/pi**2*tplmev**3*glagnem
       ndensep = 1._dl/pi**2*tplmev**3*glagnep

       if (loop.eq.1) then
         nepdiffmsw = ndensem - ndensep
       end if

       if (phie.le.17._dl) then
         thm(1) = glag1
         thm(2) = glag2
         thm(3) = glag3
         thm(4) = glag4
         thm(5) = glag5
         thm(6) = glag6
         thm(7) = glag7
         thm(11) = glag11
         thm(12) = glag12
         if (thm(12).ne.0._dl) thm(12) = 1._dl/thm(12)
       else
         thm(1) = pi**2/15._dl*tplmev**4*mev4tocgs
         thm(2) = 4._dl*thm(1)/tpl
         thm(3) = thm(1)/3._dl
         thm(4) = 0._dl
         thm(5) = 0._dl
         thm(6) = 0._dl
         thm(7) = 0._dl
         thm(11) = 0._dl
         thm(12) = 0._dl
       end if

       thm(10) = thm(1) + thm(4) + thm(8) + thm(9) &
               + rhor
       !(ref 10).
       !z = sqrt(xmelec**2/tplmev**2 + dmeeps2)
       !call bessel_eval(z)
       !thm(11) = -(z**3/tpl)*(sinh1*(3._dl*blz(1) - z*bmz(1)) - sinh2*(3._dl*blz(2) &
       !          - 2._dl*z*bmz(2)) + sinh3*(3._dl*blz(3) - 3._dl*z*bmz(3)) - sinh4 &
       !          *(3._dl*blz(4) - 4._dl*z*bmz(4)) + sinh5*(3._dl*blz(5) - 5._dl*z*bmz(5)))   !(ref 11).
       !thm(12) = z**3*(cosh1*blz(1) - 2._dl*cosh2*blz(2) &
       !          + 3._dl*cosh3*blz(3) - 4._dl*cosh4*blz(4) + 5._dl*cosh5*blz(5))   !(ref 12). 
       !write (95,*) thm(11), glag11, thm(12), glag12


       !if (rank.eq.0) write (93,*) tplmev, mt**2, dmeeps2

       !write (95,*) thm(1)/mev4tocgs/tplmev**4, (thm(1)/mev4tocgs/tplmev**4 - pi**2/15._dl)/(pi**2/15._dl)
       !write (95,*) thm(2)/mev3tocgs/tplmev**3, (thm(2)/mev3tocgs/tplmev**3 - 4._dl*pi**2/15._dl)/(4._dl*pi**2/15._dl)
       !write (95,*) thm(3)/mev4tocgs/tplmev**4, (thm(3)/mev4tocgs/tplmev**4 - pi**2/45._dl)/(pi**2/45._dl)


       return


!----------references and notes----------------------------------------
!     1)  thm(1)  = rho photon
!         (wagoner, r.v., fowler, w.a., and hoyle, f. 1967, ap. j. 148,
!          page 43, equation a2.)
!     2)  thm(2)  = d(rho photon)/d(tpl)
!     3)  thm(3)  = (p photon)/c**2
!         (wagoner, r.v., fowler, w.a., and hoyle, f. 1967,
!          page 43, equation a3.)
!     4)  thm(4)  = rho electron+positron
!         (fowler, w.a. and hoyle, f., 1964, ap. j. suppl. no. 91, 9,
!          page 281, equation b44.)
!     5)  thm(5)  = d(rho electron+positron)/d(tpl)
!     6)  thm(6)  = d(rho electron+positron)/d(phi e)
!     7)  thm(7)  = (p electron+positron)/c**2
!         (fowler, w.a. and hoyle, f., 1964, ap. j. suppl. no. 91, 9,
!          page 279, equation b27.)
!     8)  thm(8)  = rho neutrino
!                 = # neutrino species x rho electron neutrino (nondegenerate)
!                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
!     9)  thm(9)  = rho baryon
!     10) thm(10) = rho total
!                 = rho photon + rho electron+positron + rho neutrino
!                              + rho baryon
!     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
!                   d(tpl) \  2  (mc**2)**3                 /
!     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
!                   d(phi e) \  2  (mc**2)**3                 /
!     13) thm(13) = rate for n->p
!     14) thm(14) = rate for p->n


       end subroutine bbn_yes_therm_v2

