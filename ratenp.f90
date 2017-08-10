
       module ratenp

!------linkages.

!      used by - [subroutine] evolve_derivs


!------remarks.

!      Contains code to do integrals for n <-> p rates


!------modules.

       use mainvar
       use gauss
       use transvar
       use ratenp_corr
       !use jamil


       implicit none


!------temperature and mass variables.

       real(dl) mnpeps !deltamnp/tcmev
       real(dl) meps !xmelec/tcmev
       real(dl) tcmevnp !comoving temperature (MeV) for ratenp module 

       real(dl), dimension(:,:,:), allocatable :: npnuoccarray

       type(bin_scheme) nplgndr !for Gauss-Legendre Integration

       real(dl), dimension(:), allocatable :: sclgndra !scaled Gauss-Legendre abscissas
       real(dl), dimension(:), allocatable :: sclgndrw !scaled Gauss-Legendre weights


       type(bin_scheme) nplagur !for Gauss-Laguerre Integration

       real(dl), dimension(:), allocatable :: sclagura !scaled Gauss-Laguerre abscissas
       real(dl), dimension(:), allocatable :: sclagurw !scaled Gauss-Laguerre weights


       contains


!--------------------------------------------------------------------

       !subroutine ratenp_calc(tcmev,tcmpl,phie,cnorm,nuenp,for,rev)
       subroutine ratenp_calc(tcmev,tcmpl,phie,cnorm,nuenp,for,rev,ip,it)

!------linkages.

!      called by - [subroutine] evolve_derivs  
!      calls     - none
                                       

!------remarks.

!      generates rate coefficients for weak n->p and p->n reactions.
!      does not use thermal nu distributions, but rather ffact


       use renorm
#ifdef prllel
       use mpi
#endif


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

       !save


!------throughput variables.

       real(dl), intent(in) :: tcmev !Comoving temperature parameter in MeV
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter
       real(dl), intent(in) :: cnorm !normalization constant
       type(nuvar), dimension(:), intent(in) :: nuenp !\nu variable for e^- flavor
       real(dl), intent(out) :: for !n -> p rate (in s^-1)
       real(dl), intent(out) :: rev !p -> n rate (in s^-1)
       integer, intent(in) :: ip
       integer, intent(in) :: it


!------local variables.

       real(dl) :: lamfor1        !nue on n
       real(dl) :: lamfor2        !n decay
       real(dl) :: lamfor3        !e^+ on n
       real(dl) :: lamrev1        !e^- on p
       real(dl) :: lamrev2        !nuebar on p
       real(dl) :: lamrev3        !inverse n decay
       integer sumind
       real(dl) dmeeps2, dmedteps
       real(dl) dmgeps2, dmgdteps
#ifdef prllel
       integer ierr
       real(dl) :: f1prllel        !nue on n
       real(dl) :: f2prllel        !n decay
       real(dl) :: f3prllel        !e^+ on n
       real(dl) :: r1prllel        !e^- on p
       real(dl) :: r2prllel        !nuebar on p
       real(dl) :: r3prllel        !inverse n decay
#endif


!------procedure.


!------epsilon values of deltamnp and xmelec.

       mnpeps = deltamnp/tcmev
       meps = xmelec/tcmev
       tcmevnp = tcmev
       !20Aug2015 EG: the renorm corrections don't work without
       !finite-nucleon mass corrections....
       !call get_eg_mass_renorm(tcmev/tcmpl,phie,glagur,dmeeps2,dmedteps,dmgeps2,dmgdteps)
       !meps = sqrt(xmelec**2/tcmev**2 + dmeeps2/tcmpl**2)
       

       lamfor1 = 0._dl
       lamfor2 = 0._dl
       lamfor3 = 0._dl
       lamrev1 = 0._dl
       lamrev2 = 0._dl
       lamrev3 = 0._dl

       sumind = 0

!------rates with \nu_e:

       !forward n(\nu_e,e^-)p
       !lamfor1 = tcmev**5*cnorm*np_laguerre(urcafc,0._dl,tcmpl,phie,nuenp(1)%occfrac)
       lamfor1 = tcmev**5*cnorm*np_laguerre( &
                 urcafc,0._dl,tcmpl,phie,nuenp(1)%occfrac,sumind)

       !reverse p(e^-,\nu_e)n
       !lamrev1 = tcmev**5*cnorm*np_laguerre(urcarc,0._dl,tcmpl,phie,nuenp(1)%occfrac)
       lamrev1 = tcmev**5*cnorm*np_laguerre( &
                 urcarc,0._dl,tcmpl,phie,nuenp(1)%occfrac,sumind)

!------rates with \bar{\nu}_e:

       !forward n(e^+,\bar{\nu}_e)p
       !lamfor2 = tcmev**5*cnorm*np_laguerre(urcafn,mnpeps+meps,tcmpl,phie,nuenp(2)%occfrac)
       lamfor2 = tcmev**5*cnorm*np_laguerre( &
                 urcafn,mnpeps+meps,tcmpl,phie,nuenp(2)%occfrac,sumind)

       !reverse p(\bar{\nu}_e,e^+)n
       !lamrev2 = tcmev**5*cnorm*np_laguerre(urcarn,mnpeps+meps,tcmpl,phie,nuenp(2)%occfrac)
       lamrev2 = tcmev**5*cnorm*np_laguerre( &
                 urcarn,mnpeps+meps,tcmpl,phie,nuenp(2)%occfrac,sumind)

       !forward n-decay
       !lamfor3 = tcmev**5*cnorm*np_legendre(ndecay,0._dl,mnpeps-meps,tcmpl,phie,nuenp(2)%occfrac)
       lamfor3 = tcmev**5*cnorm*np_legendre( &
                 ndecay,0._dl,mnpeps-meps,tcmpl,phie,nuenp(2)%occfrac,sumind)

       !reverse inverse n-decay
       !lamrev3 = tcmev**5*cnorm*np_legendre(indecay,0._dl,mnpeps-meps,tcmpl,phie,nuenp(2)%occfrac)
       lamrev3 = tcmev**5*cnorm*np_legendre( &
                 indecay,0._dl,mnpeps-meps,tcmpl,phie,nuenp(2)%occfrac,sumind)

#ifdef prllel

       !call mpi_barrier(mpi_comm_world,ierr)

       f1prllel = 0._dl
       f2prllel = 0._dl
       f3prllel = 0._dl
       r1prllel = 0._dl
       r2prllel = 0._dl
       r3prllel = 0._dl

       call mpi_allreduce(lamfor1,f1prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(lamfor2,f2prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(lamfor3,f3prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(lamrev1,r1prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(lamrev2,r2prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(lamrev3,r3prllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       lamfor1 = f1prllel
       lamfor2 = f2prllel
       lamfor3 = f3prllel
       lamrev1 = r1prllel
       lamrev2 = r2prllel
       lamrev3 = r3prllel

#endif

!------Summed rates:
       
       for = lamfor1 + lamfor2 + lamfor3        !no unit conversions necessary
       rev = lamrev1 + lamrev2 + lamrev3


       !write (51,*) tcmev, tcmev/tcmpl, for, rev
!#ifdef prllel
!       if (rank.eq.0) then
!#endif
       !write (51,*) tcmev, lamfor1, lamfor2, lamfor3, lamrev1, lamrev2, lamrev3
       !call flush(51)
!#ifdef prllel
!       end if !(rank.eq.0)
!#endif

       return


       end subroutine ratenp_calc


!--------------------------------------------------------------------

       !real(dl) function np_legendre(func,x1,x2,tcmpl,phie,nuoccfrac)
       real(dl) function np_legendre(func,x1,x2,tcmpl,phie,nuoccfrac,sumind)

!------linkages.

!      called by - [subroutine] ratenp_calc
!      calls     - [subroutine] cscgqf


!------remarks.

!      Function to compute rate using Gauss-Legendre integration.


       implicit none


!------throughput variables.

       real(dl), external :: func !function for specific rate
       real(dl), intent(in) :: x1 !lower bound of integration
       real(dl), intent(in) :: x2 !lower bound of integration
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter
       real(dl), dimension(:), intent(in) :: nuoccfrac !\nu occupation fractions
       integer, intent(inout) :: sumind


!------local variables.

       integer i !index
       logical successflag !Flag for cscgqf
       real(dl) :: sumint !integration numbers
       real(dl)  x, val !integration numbers
       real(dl), dimension(1,1) :: tempnuoccprob !nu occupation probability at x
       real(dl) nuoccprob !nu occupation probability at x


!------procedure.

       sumint = 0._dl

       sclgndra = nplgndr%abscissas
       sclgndrw = nplgndr%weights
       call cscgqf(1,sclgndra,sclgndrw,x1,x2,successflag)

       npnuoccarray(1,1,:) = nuoccfrac(:)

       if (successflag) then
         do i=1,size(sclgndra)
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
           x = sclgndra(i)
           call interp_occfrac(eps_bins%abscissas,log(npnuoccarray),x,tempnuoccprob)
           nuoccprob = exp(tempnuoccprob(1,1))
           !nuoccprob = 1._dl/(exp(x) + 1._dl)
           val = func(x,nuoccprob,tcmpl,phie)
           sumint = sumint + val*sclgndrw(i)
#ifdef prllel
           end if
#endif
           sumind = sumind + 1
         end do !i
       end if !successflag


       np_legendre = sumint


       return


       end function np_legendre


!--------------------------------------------------------------------

       !real(dl) function np_laguerre(func,x1,tcmpl,phie,nuoccfrac)
       real(dl) function np_laguerre(func,x1,tcmpl,phie,nuoccfrac,sumind)

!------linkages.

!      called by - [subroutine] ratenp_calc
!      calls     - [subroutine] cscgqf, interp_occfrac


!------remarks.

!      Function to compute rate using Gauss-Laguerre integration.


       implicit none


!------throughput variables.

       real(dl), external :: func !function for specific rate
       real(dl), intent(in) :: x1 !lower bound of integration
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter
       real(dl), dimension(:), intent(in) :: nuoccfrac !\nu occupation fractions
       integer, intent(inout) :: sumind


!------local variables.

       integer i !index
       logical successflag !Flag for cscgqf
       real(dl) :: sumint !integration numbers
       real(dl)  x, val !integration numbers
       real(dl), dimension(1,1) :: tempnuoccprob !nu occupation probability at x
       real(dl) nuoccprob !nu occupation probability at x


!------procedure.

       sumint = 0._dl

       sclagura = nplagur%abscissas
       sclagurw = nplagur%weights
       if (x1.ne.0._dl) then
         call cscgqf(5,sclagura,sclagurw,x1,1._dl,successflag,0._dl)
       else
         successflag = .true.
       end if !x1
       sclagurw = exp(nplagur%abscissas)*sclagurw

       npnuoccarray(1,1,:) = nuoccfrac(:)

       if (successflag) then
         do i=1,size(sclagura)
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
           x = sclagura(i)
           call interp_occfrac(eps_bins%abscissas,log(npnuoccarray),x,tempnuoccprob)
           nuoccprob = exp(tempnuoccprob(1,1))
           !nuoccprob = 1._dl/(exp(x) + 1._dl)
           val = func(x,nuoccprob,tcmpl,phie)
           sumint = sumint + val*sclagurw(i)
           !if (i.eq.1) write(*,*) x, nuoccprob
#ifdef prllel
           end if
#endif
           sumind = sumind + 1
         end do !i
       end if !successflag


       np_laguerre = sumint


       return


       end function np_laguerre


!--------------------------------------------------------------------

       real(dl) function urcafc(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_laguerre
!      calls     - [function] get_coul_corr, get_zero_rad_corr


!------remarks.

!      Function to compute n(\nu_e,e^-)p integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) coulcorr !Coulomb correction
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------

!------procedure.
       
       ee = x + mnpeps
       coulcorr = get_coul_corr(-1._dl,sqrt(ee**2 - meps**2)/ee &
                               ,sqrt(ee**2 - meps**2)*tcmevnp)
       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       urcafc = coulcorr*zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                    *nuoccprob*(1._dl - 1._dl/(exp(ee*tcmpl - phie) + 1._dl)) &
                    *cvc

       


       return


       end function urcafc


!--------------------------------------------------------------------

       real(dl) function urcarc(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_laguerre
!      calls     - [function] get_coul_corr, get_zero_rad_corr


!------remarks.

!      Function to compute p(e^-,\nu_e)n integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) coulcorr !Coulomb correction
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------


!------procedure.
       
       ee = x + mnpeps
       coulcorr = get_coul_corr(-1._dl,sqrt(ee**2 - meps**2)/ee &
                               ,sqrt(ee**2 - meps**2)*tcmevnp)
       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       urcarc = coulcorr*zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                      *(1._dl - nuoccprob)/(exp(ee*tcmpl - phie) + 1._dl) &
                      *cvc


       return


       end function urcarc


!--------------------------------------------------------------------

       real(dl) function urcafn(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_laguerre
!      calls     - none


!------remarks.

!      Function to compute n(e^+,\bar{\nu}_e)p integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------


!------procedure.

       ee = x - mnpeps
       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       urcafn = zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                *(1._dl - nuoccprob)/(exp(ee*tcmpl + phie) + 1._dl) &
                *cvc


       return


       end function urcafn


!--------------------------------------------------------------------

       real(dl) function urcarn(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_laguerre
!      calls     - none


!------remarks.

!      Function to compute p(\bar{\nu}_e,e^+)n integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------


!------procedure.

       ee = x - mnpeps

       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       urcarn = zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                      *nuoccprob*(1._dl - 1._dl/(exp(ee*tcmpl + phie) + 1._dl)) &
                      *cvc


       return


       end function urcarn


!--------------------------------------------------------------------

       real(dl) function ndecay(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_legendre
!      calls     - [function] get_coul_corr, get_zero_rad_corr


!------remarks.

!      Function to compute n-decay integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) coulcorr !Coulomb correction
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------


!------procedure.
       
       ee = mnpeps - x
       coulcorr = get_coul_corr(-1._dl,sqrt(ee**2 - meps**2)/ee &
                               ,sqrt(ee**2 - meps**2)*tcmevnp)
       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       ndecay = coulcorr*zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                   *(1._dl - nuoccprob)*(1._dl - 1._dl/(exp(ee*tcmpl - phie) + 1._dl)) &
                   *cvc


       return


       end function ndecay


!--------------------------------------------------------------------

       real(dl) function indecay(x, nuoccprob, tcmpl, phie)
       

!------linkages.

!      called by - [function] np_legendre
!      calls     - [function] get_coul_corr, get_zero_rad_corr


!------remarks.

!      Function to compute inverse n-decay integrand.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: nuoccprob !\nu occupation probability
       real(dl), intent(in) :: tcmpl !Ratio of tcm/tpl
       real(dl), intent(in) :: phie !electron deg. parameter


!------local variables.

       real(dl) ee !lepton energy
       real(dl) coulcorr !Coulomb correction
       real(dl) zerorcorr !zero temperature radiative correction


!-----------------------------------------
          real(dl) cvc

          cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
          (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
          (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
          (1.43109 * (10.**(-16)) * x**6)
!-----------------------------------------


!------procedure.
       
       ee = mnpeps - x
       coulcorr = get_coul_corr(-1._dl,sqrt(ee**2 - meps**2)/ee &
                               ,sqrt(ee**2 - meps**2)*tcmevnp)
       zerorcorr = get_zero_rad_corr(sqrt(ee**2 - meps**2)/ee &
                                    , x/meps, ee/meps)

       indecay = coulcorr*zerorcorr*x**2*ee*sqrt(ee**2-meps**2) &
                   *nuoccprob/(exp(ee*tcmpl - phie) + 1._dl) &
                   *cvc


       return


       end function indecay


       end module ratenp
