
       subroutine trans_init


!------linkages.

!      called by - [subroutine] init
!      calls     - [subroutine] cgqf


!------remarks.

!      Subroutine to initialize transport variables.


!------modules.

       use transvar
       use gauss
       use trans_nunu
       use trans_nunubar
       use trans_nuer1
       use trans_nuer2
       use trans_nuepma
#ifdef prllel
       use mpi
#endif


       implicit none



!------local variables.

       integer i,n,m !indicies

       real(dl) maxepslin !Maximum epsilon value for linearly-spaced bins.
       integer lgndrbins !number of abscissas in Gauss-Legendre Integration.
       integer lagurbins !number of abscissas in Gauss-Laguerre Integration.

       logical successflag !flag for Gauss quadrature initializations.

       real(dl) parfact

       integer brnum

       integer startind !startind index for occfrac.
       real(dl) denom !denominator.

       type(nuvar), dimension(:,:), allocatable :: nutrans !FD neutrino occupation probs.

#ifdef prllel

       integer ierr !mpi variable for error instance

#endif

!------procedure.


!------set values from main_params.ini.

       open (unit=40, file='trans_params.ini', status='unknown')

       read (40,*) maxepslin
       read (40,*) lgndrbins
       read (40,*) lagurbins
       read (40,*) nuratiotol

       close (unit=40)


!------allocate arrays:

       allocate(eps_bins%abscissas(nbins)) !assign dimension of eps_bins%abscissas
       allocate(eps_bins%weights(nbins)) !assign dimension of eps_bins%weights

       allocate(glgndr%abscissas(lgndrbins)) !assign dimension of glgndr%abscissas
       allocate(glgndr%weights(lgndrbins)) !assign dimension of glgndr%weights

       allocate(glagur%abscissas(lagurbins)) !assign dimension of glagur%abscissas
       allocate(glagur%weights(lagurbins)) !assign dimension of glagur%weights


       allocate(glegai(lgndrbins)) !assign dimension of glegai
       allocate(glegwi(lgndrbins)) !assign dimension of glegwi

       allocate(glegai2(lgndrbins)) !assign dimension of glegai2
       allocate(glegwi2(lgndrbins)) !assign dimension of glegwi2
       
       allocate(glagai2(lagurbins)) !assign dimension of glagai2
       allocate(glagwi2(lagurbins)) !assign dimension of glagwi2

       allocate(glegai3(lgndrbins)) !assign dimension of glegai3
       allocate(glegwi3(lgndrbins)) !assign dimension of glegwi3

       allocate(e2a(2*lgndrbins+lagurbins)) !assign dimension of e2a
       allocate(e2w(2*lgndrbins+lagurbins)) !assign dimension of e2w

       allocate(e3a(lgndrbins)) !assign dimension of e3a
       allocate(e3w(lgndrbins)) !assign dimension of e3w

       allocate(glagai3(lagurbins)) !assign dimension of glagai3
       allocate(glagwi3(lagurbins)) !assign dimension of glagwi3

       allocate(ggenai2(max(lgndrbins,lagurbins))) !assign dimension of ggenai2
       allocate(ggenwi2(max(lgndrbins,lagurbins))) !assign dimension of ggenwi2

       allocate(e3r2a(3*lgndrbins+lagurbins)) !assign dimension of e3r2a
       allocate(e3r2w(3*lgndrbins+lagurbins)) !assign dimension of e3r2w

       allocate(algndrout(lgndrbins)) !assign dimension of algndrout
       allocate(wlgndrout(lgndrbins)) !assign dimension of wlgndrout
       
       allocate(alagurout(lagurbins)) !assign dimension of alagurout
       allocate(wlagurout(lagurbins)) !assign dimension of wlagurout

       allocate(algndrin(lgndrbins)) !assign dimension of algndrin
       allocate(wlgndrin(lgndrbins)) !assign dimension of wlgndrin

       allocate(alagurin(lagurbins)) !assign dimension of alagurin
       allocate(wlagurin(lagurbins)) !assign dimension of wlagurin

       allocate(agenout(max(lgndrbins,lagurbins))) !assign dimension of agenout
       allocate(wgenout(max(lgndrbins,lagurbins))) !assign dimension of wgenout

       allocate(agenin(max(lgndrbins,lagurbins))) !assign dimension of agenin
       allocate(wgenin(max(lgndrbins,lagurbins))) !assign dimension of wgenin

       allocate(eouta(3*lgndrbins+lagurbins)) !assign dimension of eouta
       allocate(eoutw(3*lgndrbins+lagurbins)) !assign dimension of eoutw

       allocate(ratescatt(4,2,3,2,nbins)) !assign dimension of ratescatt
       allocate(frsratescatt(4,2,3,2,nbins)) !assign dimension of frsratescatt
       allocate(rateannih(4,3,2,nbins)) !assign dimension of rateannih
       allocate(frsrateannih(4,3,2,nbins)) !assign dimension of frsrateannih

       allocate(scatt(4,2,3,2,nbins)) !assign dimension of scatt
       allocate(frsscatt(4,2,3,2,nbins)) !assign dimension of frsscatt
       allocate(annih(4,3,2,nbins)) !assign dimension of annih
       allocate(frsannih(4,3,2,nbins)) !assign dimension of frsannih

       allocate(woccfrac(3,2,nbins)) !assign dimension of woccfrac
       allocate(prec_trans(3,2,nbins)) !assign dimension of prec_trans

       allocate(nutrans(3,2)) !assign dimension of nutrans

#ifdef prllel

       allocate(rsprllel(4,2,3,2,nbins)) !assign dimension of rsprllel
       allocate(frsrsprllel(4,2,3,2,nbins)) !assign dimension of frsrsprllel
       allocate(raprllel(4,3,2,nbins)) !assign dimension of raprllel
       allocate(frsraprllel(4,3,2,nbins)) !assign dimension of frsraprllel
       allocate(dfdtprllel(3,2,nbins)) !assign dimension of dfdtprllel

#endif

!------initialize arrays:


       call cgqf(1,glgndr%abscissas,glgndr%weights,successflag)
       call cgqf(5,glagur%abscissas,glagur%weights,successflag,0._dl)


       !only works for linear spacing of bins, and nbins = 4n + 1
       !this puts a bin at zero:
       do i=1,nbins
         eps_bins%abscissas(i) = real(i-1, kind=dl)/real(nbins-1, kind=dl)*maxepslin
         if (i.eq.1) then
           eps_bins%weights(i) = 1._dl/real(nbins-1, kind=dl)*maxepslin*brends(1)
         else if (i.eq.nbins) then
           eps_bins%weights(i) = 1._dl/real(nbins-1, kind=dl)*maxepslin*brends(2)
         else
           brnum = mod(i-2,4) + 1
           eps_bins%weights(i) = 1._dl/real(nbins-1, kind=dl)*maxepslin*brcoeffs(brnum)
         end if
       end do

       maxepslinsave = maxepslin !needed for msw_epsind

       !call cgqf(5,eps_bins%abscissas,eps_bins%weights,successflag,0._dl)
       !eps_bins%weights = exp(eps_bins%abscissas)*eps_bins%weights

       nutrans(:,:)%dilfact = 1._dl
       nutrans(:,:)%mass = 0._dl
       do n=1,2
         do m=1,3
           allocate(nutrans(m,n)%occfrac(nbins))
           nutrans(m,n)%occfrac(:) = 1._dl/(exp(eps_bins%abscissas(:)) + 1._dl)
         end do !m loop
       end do !n loop

       !write (*,*) sum(nutrans(1,1)%occfrac*eps_bins%weights &
       !            *eps_bins%abscissas**3/2._dl/pi**2), &
       !            7._dl/8._dl*pi**2/30._dl



       !!02May2015 EG: calculate trans_coeffs in the future with a better algorithm.
       !!issue 1: coeffs array so large that it doesn't fit into L cache
       !!issue 2: most of the computations/subroutine calls are in the inner-most loop


       !21Apr2015 EG: calculate precision ratios.

       if (transflag) then

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
                            ,scatt,frsscatt,0._dl,1._dl,0._dl)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt

       call scattering_nuer2(nutrans,glgndr,glagur,glgndr,glagur &
                            ,scatt,frsscatt,0._dl,1._dl,0._dl)
       ratescatt = ratescatt + scatt
       frsratescatt = frsratescatt + frsscatt

       call annihilation_nuepma(nutrans,glgndr,glagur,glgndr,glagur &
                               ,annih,frsannih,0._dl,1._dl,0._dl)
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

       !prec_trans = 1.e-13_dl !for testing purposes only

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       prec_trans = 0._dl
       do i=startind,nbins
         do n=1,2
           do m=1,3
             denom = sum(frsratescatt(:,:,m,n,i)) + sum(frsrateannih(:,m,n,i))
             if (denom.gt.0._dl) then
               prec_trans(m,n,i) = abs(sum(ratescatt(:,:,m,n,i)) + sum(rateannih(:,m,n,i))) &
                                 /denom
             else
               prec_trans(m,n,i) = 1.e-13_dl
             end if
           end do
         end do
       end do


       end if !transflag


       do n=1,2
         do m=1,3
           nullify(nutrans(m,n)%occfrac)
         end do
       end do !n loop
       deallocate(nutrans)


       return


       end subroutine trans_init
