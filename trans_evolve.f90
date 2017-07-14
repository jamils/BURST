
       subroutine trans_evolve(nutrans, doccfracdt, drhonudt, tcmev, tplmev, phie)


!------linkages.

!      called by - [subroutine] bbn_evolve
!      calls     - [subroutine] cgqf


!------remarks.

!      Subroutine to calculate nu transport rates.


!------modules.

       use transvar
       use trans_nunu
       use trans_nunubar
       use trans_nuer1
       use trans_nuer2
       use trans_nuepma
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

       save


!------throughput variables.

       type(nuvar), dimension(:,:), intent(in) :: nutrans !\nu variable
       real(dl), dimension(:,:,:), intent(out) :: doccfracdt !\nu occfrac derivatives
       real(dl), intent(out) :: drhonudt !change in \nu energy density
       real(dl), intent(in) :: tcmev !comoving temp. in MeV
       real(dl), intent(in) :: tplmev !Plasma temp. in MeV
       real(dl), intent(in) :: phie !e^- deg. parameter


!------local variables.

       integer i,n,m !index
       integer startind

       real(dl) meps, tcmpl
       real(dl) numer, denom !numerator and denominator
       real(dl) test_quantity, test_ratio !testing variables
       real(dl) x, dx
       real(dl) sumrule1, sumrule2
       real(dl) sumrule3, sumrule4, sumrule5
       real(dl) sumrule1net, sumrule2net
       real(dl) sumrule1frs, sumrule2frs
       real(dl) sumrule3net, sumrule4net, sumrule5net
       real(dl) sumrule3frs, sumrule4frs, sumrule5frs
       real(dl) sumrule5n1, sumrule5n2
       real(dl) parfact
       real(dl) dmeeps2, dmedteps
       real(dl) dmgeps2, dmgdteps
       integer sumind

#ifdef prllel

       integer ierr !mpi variable for error instance
       real(dl) rhoprllel
       real(dl) netsum1p, frssum1p
       real(dl) netsum2p, frssum2p
       real(dl) netsum3p, frssum3p
       real(dl) netsum4p, frssum4p
       real(dl) netsum5p, frssum5p
       real(dl) sumrule5n1p, sumrule5n2p

#endif

!------procedure.


       meps = xmelec/tcmev
       !call get_eg_mass_renorm(tplmev,phie,glagur,dmeeps2,dmedteps,dmgeps2,dmgdteps)
       !meps = sqrt(xmelec**2 + dmeeps2*tplmev**2)/tcmev
       tcmpl = tcmev/tplmev

       !rate units: 1/(MeV^5 * s)

       ratescatt = 0._dl
       frsratescatt = 0._dl
       rateannih = 0._dl
       frsrateannih = 0._dl

       call scattering_nunu(nutrans,glgndr,scatt,frsscatt)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt

       call scattering_nunubar(nutrans,glgndr,scatt,frsscatt,annih,frsannih)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt
       rateannih = rateannih + annih
       frsrateannih = frsrateannih + frsannih

       call scattering_nuer1(nutrans,glgndr,glagur,glgndr &
                            ,scatt,frsscatt,meps,tcmpl,phie)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt

       call scattering_nuer2(nutrans,glgndr,glagur,glgndr,glagur &
                            ,scatt,frsscatt,meps,tcmpl,phie)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt

       call annihilation_nuepma(nutrans,glgndr,glagur,glgndr,glagur &
                               ,annih,frsannih,meps,tcmpl,phie)
       rateannih = rateannih + annih
       frsrateannih = frsrateannih + frsannih

#ifdef prllel

       call mpi_barrier(mpi_comm_world,ierr)

       rsprllel = 0._dl
       frsrsprllel = 0._dl
       raprllel = 0._dl
       frsraprllel = 0._dl

       call mpi_allreduce(ratescatt,rsprllel,4*2*3*2*nbins, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(frsratescatt,frsrsprllel,4*2*3*2*nbins, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(rateannih,raprllel,4*3*2*nbins, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(frsrateannih,frsraprllel,4*3*2*nbins, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       ratescatt = rsprllel
       frsratescatt = frsrsprllel
       rateannih = raprllel
       frsrateannih = frsraprllel

#endif


       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       doccfracdt = 0._dl
       drhonudt = 0._dl

       sumrule1 = 0._dl
       sumrule2 = 0._dl
       sumrule1net = 0._dl
       sumrule1frs = 0._dl
       sumrule2net = 0._dl
       sumrule2frs = 0._dl
       sumrule3net = 0._dl
       sumrule3frs = 0._dl
       sumrule4net = 0._dl
       sumrule4frs = 0._dl
       sumrule5net = 0._dl
       sumrule5frs = 0._dl
       sumrule5n1 = 0._dl
       sumrule5n2 = 0._dl

       sumind = 0

       do i=startind,nbins
         x = eps_bins%abscissas(i)
         dx = eps_bins%weights(i)
         parfact = 1._dl
         do n=1,2
           do m=1,3
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             numer = sum(ratescatt(:,:,m,n,i)) + sum(rateannih(:,m,n,i))
             denom = sum(frsratescatt(:,:,m,n,i)) + sum(frsrateannih(:,m,n,i))
             if (denom.ne.0._dl) then
               test_quantity = abs(numer)/denom
             else
               test_quantity = 0._dl
             end if
             test_ratio = test_quantity/prec_trans(m,n,i)
             !conditioning of rates:
             if (test_ratio.ge.nuratiotol) then !implemented from trans_params.ini
               doccfracdt(m,n,i) = tcmev**5*numer !units of s^{-1}
               !Change in plasma energy
               drhonudt = drhonudt + 1._dl/2._dl/pi**2*tcmev**9*dx*x**3* &
                          (sum(ratescatt(4,:,m,n,i)) + rateannih(4,m,n,i)) !units of MeV^4/s

               !Sum rule tests for error monitoring:
               sumrule1net = sumrule1net + 1._dl/2._dl/pi**2*dx*x**2* &
                          (sum(ratescatt(1:3,:,m,n,i)) + sum(rateannih(1:3,m,n,i)))
               sumrule1frs = sumrule1frs + 1._dl/2._dl/pi**2*dx*x**2* &
                          (sum(frsratescatt(1:3,:,m,n,i)) + sum(frsrateannih(1:3,m,n,i)))
               sumrule2net = sumrule2net + 1._dl/2._dl/pi**2*dx*x**3* &
                          (sum(ratescatt(1:3,:,m,n,i)) + sum(rateannih(1:3,m,n,i)))
               sumrule2frs = sumrule2frs + 1._dl/2._dl/pi**2*dx*x**3* &
                          (sum(frsratescatt(1:3,:,m,n,i)) + sum(frsrateannih(1:3,m,n,i)))
               sumrule3net = sumrule3net + 1._dl/2._dl/pi**2*dx*x**2* &
                          parfact*sum(rateannih(1:4,m,n,i))
               sumrule3frs = sumrule3frs + 1._dl/2._dl/pi**2*dx*x**2* &
                          sum(frsrateannih(1:4,m,n,i))
               sumrule4net = sumrule4net + 1._dl/2._dl/pi**2*dx*x**2* &
                          parfact*sum(rateannih(1:3,m,n,i))
               sumrule4frs = sumrule4frs + 1._dl/2._dl/pi**2*dx*x**2* &
                          sum(frsrateannih(1:3,m,n,i))
               sumrule5net = sumrule5net + 1._dl/2._dl/pi**2*dx*x**2* &
                          parfact*rateannih(4,m,n,i)
               sumrule5frs = sumrule5frs + 1._dl/2._dl/pi**2*dx*x**2* &
                          frsrateannih(4,m,n,i)
               if (n.eq.1) then
                 sumrule5n1 = sumrule5n1 + 1._dl/2._dl/pi**2*dx*x**2* &
                            rateannih(4,m,n,i)
               else
                 sumrule5n2 = sumrule5n2 + 1._dl/2._dl/pi**2*dx*x**2* &
                            rateannih(4,m,n,i)
               end if
             end if
#ifdef prllel
           end if
#endif
             sumind = sumind + 1
           end do !m loop
           parfact = -1._dl
         end do !n loop
       end do !i loop

#ifdef prllel

       call mpi_barrier(mpi_comm_world,ierr)

       rhoprllel = 0._dl
       dfdtprllel = 0._dl

       netsum1p = 0._dl
       frssum1p = 0._dl
       netsum2p = 0._dl
       frssum2p = 0._dl
       netsum3p = 0._dl
       frssum3p = 0._dl
       netsum4p = 0._dl
       frssum4p = 0._dl
       netsum5p = 0._dl
       frssum5p = 0._dl
       sumrule5n1p = 0._dl
       sumrule5n2p = 0._dl

       call mpi_allreduce(doccfracdt,dfdtprllel,3*2*nbins, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(drhonudt,rhoprllel,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       call mpi_allreduce(sumrule1net,netsum1p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule1frs,frssum1p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule2net,netsum2p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule2frs,frssum2p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule3net,netsum3p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule3frs,frssum3p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule4net,netsum4p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule4frs,frssum4p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule5net,netsum5p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule5frs,frssum5p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule5n1,sumrule5n1p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(sumrule5n2,sumrule5n2p,1, &
            mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       doccfracdt = dfdtprllel
       drhonudt = rhoprllel

       sumrule1net = netsum1p
       sumrule1frs = frssum1p
       sumrule2net = netsum2p
       sumrule2frs = frssum2p
       sumrule3net = netsum3p
       sumrule3frs = frssum3p
       sumrule4net = netsum4p
       sumrule4frs = frssum4p
       sumrule5net = netsum5p
       sumrule5frs = frssum5p
       sumrule5n1 = sumrule5n1p
       sumrule5n2 = sumrule5n2p

#endif

       sumrule1 = sumrule1net/sumrule1frs
       sumrule2 = sumrule2net/sumrule2frs
       sumrule3 = sumrule3net/sumrule3frs
       sumrule4 = sumrule4net/sumrule4frs
       sumrule5 = sumrule5net/sumrule5frs

#ifdef prllel
       if (rank.eq.mod(0,nprocs)) then
#endif
       write (52,*) tcmev, tplmev, sumrule1, sumrule2, sumrule1net, sumrule1frs, sumrule2net, sumrule2frs
       call flush(52)
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.mod(1,nprocs)) then
#endif
       !write (53,*) tcmev, sumrule3, sumrule4, sumrule5 
       !write (53,*) tcmev, sumrule3, sumrule3net, sumrule3frs 
       !write (53,*) tcmev, &
       !sumrule5n1,sumrule5n2,sumrule5net,sumrule5frs,sumrule5
       !write (53,*) tcmev, frsratescatt(4,1,1,1,:)
       !call flush(53)
#ifdef prllel
       end if
#endif


       return


       end subroutine trans_evolve
