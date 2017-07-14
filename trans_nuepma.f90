
       module trans_nuepma


!------linkages.

!      used by main

!------remarks.

!      Module for electron positron annihilation into neutrinos.


       use mainvar
       use transvar
       use interp


       implicit none


       contains


!-----------------------------------------------------------

       !subroutine annihilation_nuepma(nutrans,eps_bins,out_bins_lgndr,out_bins_lagur &
       !                            ,in_bins_lgndr,in_bins_lagur,net,frs,meps,tcmpl,phie)
       subroutine annihilation_nuepma(nutrans,out_bins_lgndr,out_bins_lagur &
                                   ,in_bins_lgndr,in_bins_lagur,net,frs,meps,tcmpl,phie)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] none
!                  [function] none


!-----remarks.

!     Calculates annihilation rates for \nu-\bar{\nu} into e^\pm.
!     out_bins_lgndr and out_bins_lagur must have the same size


!------modules.

       use gauss
#ifdef prllel
       use mpi
#endif


       implicit none


       save


!------throughput variables.

       type (nuvar), intent(in), dimension(:,:) :: nutrans !nu occupation fractions and info.
       !type (bin_scheme), intent(in) :: eps_bins !epsilon abscissas and weights
       type (bin_scheme), intent(in) :: out_bins_lgndr !Gauss-Legendre abscissas and weights for outer integral
       type (bin_scheme), intent(in) :: out_bins_lagur !Gauss-Laguerre abscissas and weights for outer integral
       type (bin_scheme), intent(in) :: in_bins_lgndr !Gauss-Legendre abscissas and weights for inner integral
       type (bin_scheme), intent(in) :: in_bins_lagur !Gauss-Laguerre abscissas and weights for inner integral
       real(dl), intent(out), dimension(:,:,:,:) :: net !net annihilation rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:) :: frs !forward-reverse summed annihilation rates/tnmev**5
       real(dl), intent(in) :: meps !m_e/T_{cm}
       real(dl), intent(in) :: tcmpl !T_{cm}/T_{pl}
       real(dl), intent(in) :: phie !electron degeneracy parameter


!------local variables.

       integer i,j,k,m,mm,n
       !integer nbins
       integer sumind

       real(dl) p1, e2, q2, e3, q3, p4, eout, qout, ein, qin
       real(dl) prefac
       real(dl) kinemfact1
       real(dl), dimension(3,2) :: kinemfactin
       real(dl), dimension(3,2) :: kinemfactinswap
       real(dl) parfact
       real(dl) pcase1, pcase2, pcase3
       real(dl) ecut1, ecut2
       real(dl) etrans1, etrans2
       real(dl) elim1, elim2
       real(dl) ediff

       real(dl), dimension(3,2) :: logwoccfrac

       real(dl), dimension(4,2,4) :: woccarray !array for calculating F

       logical successflag, ecutflag

       !real(dl), dimension(:), allocatable :: algndrout
       !real(dl), dimension(:), allocatable :: wlgndrout

       !real(dl), dimension(:), allocatable :: alagurout
       !real(dl), dimension(:), allocatable :: wlagurout

       !real(dl), dimension(:), allocatable :: algndrin
       !real(dl), dimension(:), allocatable :: wlgndrin

       !real(dl), dimension(:), allocatable :: alagurin
       !real(dl), dimension(:), allocatable :: wlagurin

       !real(dl), dimension(:), allocatable :: agenout
       !real(dl), dimension(:), allocatable :: wgenout

       !real(dl), dimension(:), allocatable :: agenin
       !real(dl), dimension(:), allocatable :: wgenin

       !real(dl), dimension(:), allocatable :: eouta
       !real(dl), dimension(:), allocatable :: eoutw

       real(dl), dimension(3,2) :: weightin
       real(dl), dimension(3,2) :: weightinswap

       real(dl), dimension(3,2) :: net1
       real(dl), dimension(3,2) :: frs1
       real(dl), dimension(3,2) :: net2
       real(dl), dimension(3,2) :: frs2
       real(dl), dimension(3,2) :: net2swap
       real(dl), dimension(3,2) :: frs2swap

       real(dl) botlim, toplim
       real(dl) tempint

       integer startind !starting index
#ifdef prllel

       integer ierr !mpi variable for error instance

#endif


!------Testing quantities.

       !real(dl), dimension(:,:), allocatable :: testsum
       !real(dl) tempint
       !real(dl) lowsum, highsum
       !real(dl) testsum2, woccfrac_exact
       !real(dl), dimension(3,100,3) :: allarray
       !integer jj


!------procedure.


       !allocate(algndrout(size(out_bins_lgndr%abscissas))) !assign dimension of algndrout
       !allocate(wlgndrout(size(out_bins_lgndr%abscissas))) !assign dimension of wlgndrout
       
       !allocate(alagurout(size(out_bins_lagur%abscissas))) !assign dimension of alagurout
       !allocate(wlagurout(size(out_bins_lagur%abscissas))) !assign dimension of wlagurout

       !allocate(algndrin(size(in_bins_lgndr%abscissas))) !assign dimension of algndrin
       !allocate(wlgndrin(size(in_bins_lgndr%abscissas))) !assign dimension of wlgndrin

       !allocate(alagurin(size(in_bins_lagur%abscissas))) !assign dimension of alagurin
       !allocate(wlagurin(size(in_bins_lagur%abscissas))) !assign dimension of wlagurin

       !allocate(agenout(max(size(out_bins_lgndr%abscissas),size(out_bins_lagur%abscissas)))) !assign dimension of agenout
       !allocate(wgenout(max(size(out_bins_lgndr%abscissas),size(out_bins_lagur%abscissas)))) !assign dimension of wgenout

       !allocate(agenin(max(size(in_bins_lgndr%abscissas),size(in_bins_lagur%abscissas)))) !assign dimension of agenin
       !allocate(wgenin(max(size(in_bins_lgndr%abscissas),size(in_bins_lagur%abscissas)))) !assign dimension of wgenin

       !allocate(eouta(3*size(algndrout)+size(alagurout))) !assign dimension of eouta
       !allocate(eoutw(3*size(algndrout)+size(alagurout))) !assign dimension of eoutw


       !nbins = size(nutrans(1,1)%occfrac)

       !allocate(testsum(nbins,3*size(out_bins_lgndr%abscissas)+size(out_bins_lagur%abscissas))) !assign dimension of testsum


       do n=1,2
         do m=1,3 !don't use sterile occupation fractions
           woccfrac(m,n,:) = nutrans(m,n)%occfrac(:)
         end do
       end do
       woccarray = 0._dl
       eouta = 0._dl
       eoutw = 0._dl


       prefac = gfhbar/16._dl/(2._dl*pi)**3
       !testsum = 0._dl
       !testsum2 = 0._dl

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       !p_1 cases:
       pcase1 = meps/2._dl
       pcase2 = meps*(1._dl + sqrt(5._dl))/4._dl
       pcase3 = meps


       net = 0._dl
       frs = 0._dl
       sumind = 1
       !write (116,'(A,1pe12.5)') 'm_e = ', meps
       !write (119,'(A,1pe12.5)') 'm_e = ', meps
       !write (120,'(A,1pe12.5)') 'm_e = ', meps
       !write (121,'(A,1pe12.5)') 'm_e = ', meps
       !write (118,'(A,1pe12.5)') 'm_e = ', meps
       do i=startind,nbins !outermost loop over p_1
       !do i=startind,nbins,10 !outermost loop over p_1
       !do i=17,17 !outermost loop over p_1
         !Energy and occupation probabilities for first particle:
         p1 = eps_bins%abscissas(i)
         woccarray(1:3,:,1) = woccfrac(1:3,:,i)
         !Outer integral cut points:
         ecut1 = p1 + meps**2/4._dl/p1 !=\ecuti, when \etransi is used
         ecut2 = 0.5_dl*(2._dl*p1 - meps + meps**2/(2._dl*p1 - meps)) !=\ecutii, when \elimi and \etransii are used

         !write (116,*)
         !write (116,'(A,1pe12.5)') 'p_1 = ', p1
         !write (119,*)
         !write (119,'(A,1pe12.5)') 'p_1 = ', p1
         !write (120,*)
         !write (120,'(A,1pe12.5)') 'p_1 = ', p1
         !write (121,*)
         !write (121,'(A,1pe12.5)') 'p_1 = ', p1
         !write (116,*) 'p_1 = ', p1, ecut1
         !write (118,*)
         !write (118,'(A,1pe12.5)') 'p1 = ', p1

         !p1/meps < 1/2:
         if (p1.lt.pcase1) then !ecut2 unneccessary

           !Make the abscissas and weights for \int_{\ecuti}^{\infty} dE_out

           !Gauss Laguerre (lagur) interval:
           alagurout = out_bins_lagur%abscissas
           wlagurout = out_bins_lagur%weights
           call cscgqf(5,alagurout,wlagurout,ecut1,1._dl,successflag,0._dl)
           !a and w for \int_{\ecuti}^{\infty} dE_out
           eouta(1:size(alagurout)) = alagurout
           eoutw(1:size(alagurout)) = exp(out_bins_lagur%abscissas)*wlagurout

           do j=1,size(alagurout)
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             !Outer integral energy and momentum:
             eout = eouta(j)
             qout = sqrt(eout**2 - meps**2)
             !Inner integral transition and limit points:
             elim1 = 0.5_dl*(2._dl*p1 - eout + qout + meps**2/(2._dl*p1 - eout + qout))
             !Occupation probabilities for 1st lepton (not neccessarily E_2):
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,2) = fd_equil_calc(eout,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             !Kinematic factor for outer integral:
             kinemfact1 = eoutw(j)/p1**2
             !Initialize temporary arrays for net and frs rates:
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'eout = ', eout
             !write (119,*)
             !write (119,*) 'eout = ', eout
             !write (120,*)
             !write (120,*) 'eout = ', eout
             !write (121,*)
             !write (121,*) 'eout = ', eout
             !write (116,*) 'eout = ', eout, elim1, j
             !tempint = 0._dl
             !write (*,*) elim2, eout, elim1

             !\int_{\ecuti}^{\infty} dE_out:

             !if statement protects against machine imprecision,
             !set \int dE_in to zero if \elimi is large negative number
             if (elim1.gt.0._dl) then 
               !a and w for \int_{\elimi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,elim1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 !Inner integral energy and momentum:
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 !Occupation probabilities for lepton B (not neccessarily E_3):
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 !Energy and occupation probabilities for fourth particle:
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 !Generic limits for \int dy:
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 !when E_out = E_2:
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 !when E_out = E_3:
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 !Calculate total net and for rates for both integrals:
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimi}^{\infty} dE_in
             end if !(elim1.gt.0._dl)


             !if (i.eq.2) write (200+rank,*) eout, tempint,net(4,1,1,2)

             !No swap variable needed as the integral is identical for either E_2 or E_3:
             net(4,:,:,i) = net(4,:,:,i) + net1(:,:)*kinemfact1
             frs(4,:,:,i) = frs(4,:,:,i) + frs1(:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !1/2 < p1/meps < (\sqrt{5} + 1)/4:
         else if ((p1.lt.pcase2).and.(p1.ge.pcase1)) then !ecut1 < ecut2

           !Make the abscissas and weights for \int_{m_e}^{\infty} dE_out

           !Gauss Legndre (lgndr) interval:
           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,meps,ecut1,successflag)
           !a and w for \int_{m_e}^{\ecuti} dE_out
           if (successflag) then
             eouta(1:size(algndrout)) = algndrout
             eoutw(1:size(algndrout)) = wlgndrout
           else
             eouta(1:size(algndrout)) = meps
             eoutw(1:size(algndrout)) = 0._dl
           end if !successflag

           !put in general abscissas and weights:
           !a and w for \int_{\ecuti}^{\ecutii} dE_out
           if ((ecut2-ecut1).lt.200._dl) then !Use Legendre
             agenout = out_bins_lgndr%abscissas
             wgenout = out_bins_lgndr%weights
             call cscgqf(1,agenout,wgenout,ecut1,ecut2,successflag)
             if (successflag) then
               eouta(1+size(algndrout):2*size(algndrout)) = agenout
               eoutw(1+size(algndrout):2*size(algndrout)) = wgenout
             else
               eouta(1+size(algndrout):2*size(algndrout)) = ecut1
               eoutw(1+size(algndrout):2*size(algndrout)) = 0._dl
             end if !successflag
           else !Use Laguerre
             agenout = out_bins_lagur%abscissas
             wgenout = out_bins_lagur%weights
             call cscgqf(5,agenout,wgenout,ecut1,1._dl,successflag,0._dl)
             wgenout = exp(out_bins_lagur%abscissas)*wgenout
             eouta(1+size(algndrout):2*size(algndrout)) = agenout
             eoutw(1+size(algndrout):2*size(algndrout)) = wgenout
           end if

           !Gauss Laguerre (lagur) interval:
           alagurout = out_bins_lagur%abscissas
           wlagurout = out_bins_lagur%weights
           call cscgqf(5,alagurout,wlagurout,ecut2,1._dl,successflag,0._dl)
           !a and w for \int_{\ecutii}^{\infty} dE_out
           eouta(1+2*size(algndrout):2*size(algndrout)+size(alagurout)) = alagurout
           eoutw(1+2*size(algndrout):2*size(algndrout)+size(alagurout)) = exp(out_bins_lagur%abscissas)*wlagurout

           do j=1,2*size(algndrout)+size(alagurout)
           !do j=10,2*size(algndrout)+size(alagurout),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             eout = eouta(j)
             qout = sqrt(eout**2 - meps**2)
             etrans1 = 0.5_dl*(2._dl*p1 - eout - qout + meps**2/(2._dl*p1 - eout - qout))
             !etrans2 = 0.5_dl*(2._dl*p1 - eout + qout + meps**2/(2._dl*p1 - eout + qout))
             etrans2 = max(meps,0.5_dl*(2._dl*p1 - eout + qout + meps**2/(2._dl*p1 - eout + qout)))
             !if (etrans2.lt.meps) then
             !  write (*,*) i,j,meps, etrans2, p1, eout
             !  stop
             !end if
             elim1 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,2) = fd_equil_calc(eout,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = eoutw(j)/p1**2
             !Initialize temporary arrays for net and frs rates:
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'eout = ', eout
             !write (119,*)
             !write (119,*) 'eout = ', eout
             !write (120,*)
             !write (120,*) 'eout = ', eout
             !write (121,*)
             !write (121,*) 'eout = ', eout
             !tempint = 0._dl
             !write (*,*) elim2, eout, elim1

             !\int_{m_e}^{\ecuti} dE_out:
             if (j.le.size(algndrout)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransi is large negative number
               ediff = etrans1 - elim1
               !put in genin:
               !a and w frs \int_{\elimi}^{\etransi} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,elim1,etrans1,successflag)
                 if (.not.successflag) then
                   write (58,*) 'case211',i,j, ecut1, ecut2, elim1, etrans1
                   agenin = elim1
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,elim1,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !if (isnan(net1(1,1))) write (58,*) i,j,k, kinemfactin(1,1), kinemfactinswap(1,1)
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimi}^{\etransi} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransi is large negative number
               if (etrans1.gt.0._dl) then
               !a and w for \int_{\etransi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 if (isnan(net1(1,1))) write (58,*) i,j,k, kinemfactin(1,1), kinemfactinswap(1,1)
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimi}^{\infty} dE_in
               end if !(etrans1.gt.0._dl)

              
             !\int_{\ecuti}^{\ecutii} dE_out:
             else if ((j.le.2*size(algndrout)).and.(j.ge.size(algndrout)+1)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \elimi is large negative number
               if (elim1.gt.0._dl) then
               !a and w for \int_{\elimi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,elim1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 if (isnan(net1(1,1))) write (58,*) i,j,k, kinemfactin(1,1), kinemfactinswap(1,1)
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimi}^{\infty} dE_in
               end if !(elim1.gt.0._dl)


             !\int_{\ecutii}^{\infty} dE_out:
             else if ((j.le.2*size(algndrout)+size(alagurout)).and.(j.ge.2*size(algndrout)+1)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransi is large negative number
               ediff = etrans2 - meps
               !put in genin:
               !a and w frs \int_{m_e}^{\etransii} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,meps,etrans2,successflag)
                 if (.not.successflag) then
                   write (58,*) 'case231',i,j, ecut1, ecut2, meps, etrans2
                   agenin = elim1
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,meps,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 !qin = sqrt(ein**2 - meps**2) !q_in momentum
                 qin = dsqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !if (isnan(net1(1,1))) then
                 !  write (90,*) woccarray(:,:,4)
                 !  write (90,*) woccfrac(1,1,:)
                 !  call flush(90)
                 !  stop
                 !end if
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransii is large negative number
               if (etrans2.gt.0._dl) then
               !a and w for \int_{\etransii}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans2,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 if (isnan(net1(1,1))) write (58,*) i,j,k, kinemfactin(1,1), kinemfactinswap(1,1)
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\infty} dE_in
               end if !(etrans2.gt.0._dl)


             end if !\int dE_out


             !write (116,'(/)')
             !write (118,*) eout, tempint

             net(4,:,:,i) = net(4,:,:,i) + net1(:,:)*kinemfact1
             frs(4,:,:,i) = frs(4,:,:,i) + frs1(:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !(\sqrt{5} + 1)/4 < p1/meps < 1:
         else if ((p1.lt.pcase3).and.(p1.ge.pcase2)) then !p1 < ecut2 < ecut1

           !Make the abscissas and weights for \int_{m_e}^{\infty} dE_out

           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,meps,ecut2,successflag)
           !a and w for \int_{m_e}^{\ecutii} dE_out
           eouta(1:size(algndrout)) = algndrout
           eoutw(1:size(algndrout)) = wlgndrout

           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,ecut2,ecut1,successflag)
           !a and w for \int_{\ecutii}^{\ecuti} dE_out
           eouta(1+size(algndrout):2*size(algndrout)) = algndrout
           eoutw(1+size(algndrout):2*size(algndrout)) = wlgndrout

           !Gauss Laguerre (lagur) interval:
           alagurout = out_bins_lagur%abscissas
           wlagurout = out_bins_lagur%weights
           call cscgqf(5,alagurout,wlagurout,ecut1,1._dl,successflag,0._dl)
           !a and w for \int_{\ecuti}^{\infty} dE_out
           eouta(1+2*size(algndrout):2*size(algndrout)+size(alagurout)) = alagurout
           eoutw(1+2*size(algndrout):2*size(algndrout)+size(alagurout)) = exp(out_bins_lagur%abscissas)*wlagurout

           do j=1,2*size(algndrout)+size(alagurout)
           !do j=10,2*size(algndrout)+size(alagurout),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             eout = eouta(j)
             qout = sqrt(eout**2 - meps**2)
             etrans1 = 0.5_dl*(2._dl*p1 - eout - qout + meps**2/(2._dl*p1 - eout - qout))
             etrans2 = 0.5_dl*(2._dl*p1 - eout + qout + meps**2/(2._dl*p1 - eout + qout))
             elim1 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,2) = fd_equil_calc(eout,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = eoutw(j)/p1**2
             !Initialize temporary arrays for net and frs rates:
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'eout = ', eout
             !write (119,*)
             !write (119,*) 'eout = ', eout
             !write (120,*)
             !write (120,*) 'eout = ', eout
             !write (121,*)
             !write (121,*) 'eout = ', eout
             !write (116,*) 'eout = ', eout, elim1, j
             !tempint = 0._dl
             !write (*,*) elim2, eout, elim1


             !\int_{m_e}^{\ecutii} dE_out:
             if (j.le.size(algndrout)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransi is large negative number
               ediff = etrans1 - elim1
               !put in genin:
               !a and w frs \int_{\elimi}^{\etransi} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,elim1,etrans1,successflag)
                 if (.not.successflag) then
                   write (58,*) i,j, ecut1, ecut2, elim1, etrans1
                   agenin = elim1
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,elim1,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimi}^{\etransi} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransi is large negative number
               if (etrans1.gt.0._dl) then
               !a and w for \int_{\etransi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransi}^{\infty} dE_in
               end if !(etrans1.gt.0._dl)


             !\int_{\ecutii}^{\ecuti} dE_out:
             else if ((j.le.2*size(algndrout)).and.(j.ge.size(algndrout)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_in
               algndrin = in_bins_lgndr%abscissas
               wlgndrin = in_bins_lgndr%weights
               call cscgqf(1,algndrin,wlgndrin,meps,etrans2,successflag)
               if (.not.successflag) then
                 !write (58,*) i,j, ecut1, ecut2, meps, etrans2
                 algndrin = meps
                 wlgndrin = 0._dl
               end if !successflag
               do k=1,size(algndrin)
                 ein = algndrin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlgndrin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlgndrin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransi is large negative number
               ediff = etrans1 - etrans2
               !put in genin:
               !a and w frs \int_{\etransii}^{\etransi} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,etrans2,etrans1,successflag)
                 if (.not.successflag) then
                   write (58,*) i,j, ecut1, ecut2, etrans1, etrans2
                   agenin = etrans2
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,etrans2,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\etransi} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransi is large negative number
               if (etrans1.gt.0._dl) then
               !a and w for \int_{\etransi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransi}^{\infty} dE_in
               end if !(etrans1.gt.0._dl)


             !\int_{\ecuti}^{\infty} dE_out:
             else if ((j.le.2*size(algndrout)+size(alagurout)).and.(j.ge.2*size(algndrout)+1)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransii is large negative number
               ediff = etrans2 - meps
               !put in genin:
               !a and w frs \int_{m_e}^{\etransii} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,meps,etrans2,successflag)
                 if (.not.successflag) then
                   !write (58,*) i,j, ecut1, ecut2, meps, etrans2
                   agenin = meps
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,meps,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransii is large negative number
               if (etrans2.gt.0._dl) then
               !a and w for \int_{\etransii}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans2,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\infty} dE_in
               end if !(etrans2.gt.0._dl)


             end if !\int dE_out


             !write (118,*) eout, tempint

             net(4,:,:,i) = net(4,:,:,i) + net1(:,:)*kinemfact1
             frs(4,:,:,i) = frs(4,:,:,i) + frs1(:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !1 < p1/meps:
         else !ecut2 < p1 < ecut1

           ecutflag = .true.

           !Make the abscissas and weights for \int_{m_e}^{\infty} dE_out

           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,meps,ecut2,successflag)
           !a and w for \int_{m_e}^{\ecutii} dE_out
           if (successflag) then
             eouta(1:size(algndrout)) = algndrout
             eoutw(1:size(algndrout)) = wlgndrout
           else !09Apr2015 EG: in case p_1 == m_e => ecut2 = m_e;
             eouta(1:size(algndrout)) = meps
             eoutw(1:size(algndrout)) = 0._dl
             ecutflag = .false.
           end if

           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,ecut2,p1,successflag)
           !a and w for \int_{\ecutii}^{p_1} dE_out
           if (successflag) then
             eouta(1+size(algndrout):2*size(algndrout)) = algndrout
             eoutw(1+size(algndrout):2*size(algndrout)) = wlgndrout
           else !09Apr2015 EG: in case p_1 == m_e => ecut2 = m_e;
             eouta(1+size(algndrout):2*size(algndrout)) = meps
             eoutw(1+size(algndrout):2*size(algndrout)) = 0._dl
             ecutflag = .false.
           end if

           algndrout = out_bins_lgndr%abscissas
           wlgndrout = out_bins_lgndr%weights
           call cscgqf(1,algndrout,wlgndrout,p1,ecut1,successflag)
           !a and w for \int_{p_1}^{\ecuti} dE_out
           if (successflag) then
             eouta(1+2*size(algndrout):3*size(algndrout)) = algndrout
             eoutw(1+2*size(algndrout):3*size(algndrout)) = wlgndrout
           else !10Apr2015 EG: in case m_e = 0 => ecut1 = p1;
             eouta(1+2*size(algndrout):3*size(algndrout)) = p1
             eoutw(1+2*size(algndrout):3*size(algndrout)) = 0._dl
             ecutflag = .false.
           end if

           !Gauss Laguerre (lagur) interval:
           alagurout = out_bins_lagur%abscissas
           wlagurout = out_bins_lagur%weights
           call cscgqf(5,alagurout,wlagurout,ecut1,1._dl,successflag,0._dl)
           !a and w for \int_{\ecuti}^{\infty} dE_out
           eouta(1+3*size(algndrout):3*size(algndrout)+size(alagurout)) = alagurout
           eoutw(1+3*size(algndrout):3*size(algndrout)+size(alagurout)) = exp(out_bins_lagur%abscissas)*wlagurout


           do j=1,3*size(algndrout)+size(alagurout)
           !do j=10,3*size(algndrout)+size(alagurout),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             eout = eouta(j)
             qout = sqrt(eout**2 - meps**2)
             etrans1 = 0.5_dl*(2._dl*p1 - eout - qout + meps**2/(2._dl*p1 - eout - qout))
             etrans2 = 0.5_dl*(2._dl*p1 - eout + qout + meps**2/(2._dl*p1 - eout + qout))
             elim1 = etrans2
             elim2 = etrans1
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,2) = fd_equil_calc(eout,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = eoutw(j)/p1**2
             !Initialize temporary arrays for net and frs rates:
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'eout = ', eout
             !write (119,*)
             !write (119,*) 'eout = ', eout
             !write (120,*)
             !write (120,*) 'eout = ', eout
             !write (121,*)
             !write (121,*) 'eout = ', eout
             !write (122,*)
             !write (122,*) etrans1,etrans2
             !write (116,*) 'eout = ', eout, elim1, j
             !tempint = 0._dl
             !write (*,*) elim2, eout, elim1


             !\int_{m_e}^{\ecutii} dE_out:
             if (j.le.size(algndrout) &
                       .and.ecutflag) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransii is large negative number
               ediff = etrans2 - elim2
               !put in genin:
               !a and w frs \int_{\elimii}^{\etransii} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,elim2,etrans2,successflag)
                 if (.not.successflag) then
                   !write (58,*) i,j, ecut1, ecut2, elim2, etrans2
                   agenin = elim2
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,elim2,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\elimii}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransii is large negative number
               if (etrans2.gt.0._dl) then
               !a and w for \int_{\etransii}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans2,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\infty} dE_in
               end if !(etrans2.gt.0._dl)


             !\int_{\ecutii}^{p_1} dE_out:
             else if ((j.le.2*size(algndrout)).and.(j.ge.size(algndrout)+1) &
                       .and.ecutflag) then

               !a and w for \int_{m_e}^{\etransi} dE_in
               algndrin = in_bins_lgndr%abscissas
               wlgndrin = in_bins_lgndr%weights
               call cscgqf(1,algndrin,wlgndrin,meps,etrans1,successflag)
               if (.not.successflag) then
                 !write (58,*) i,j, ecut1, ecut2, meps, etrans1
                 algndrin = meps
                 wlgndrin = 0._dl
               end if !successflag
               do k=1,size(algndrin)
                 ein = algndrin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlgndrin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlgndrin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransi} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransii is large negative number
               ediff = etrans2 - etrans1
               !put in genin:
               !a and w frs \int_{\etransi}^{\etransii} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,etrans1,etrans2,successflag)
                 if (.not.successflag) then
                   !write (58,*) i,j, ecut1, ecut2, etrans1, etrans2
                   agenin = etrans1
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,etrans1,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransi}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransii is large negative number
               if (etrans2.gt.0._dl) then
               !a and w for \int_{\etransii}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans2,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\infty} dE_in
               end if !(etrans2.gt.0._dl)


             !\int_{p_1}^{\ecuti} dE_out:
             else if ((j.le.3*size(algndrout)).and.(j.ge.2*size(algndrout)+1) &
                       .and.ecutflag) then

               !a and w for \int_{m_e}^{\etransii} dE_in
               algndrin = in_bins_lgndr%abscissas
               wlgndrin = in_bins_lgndr%weights
               call cscgqf(1,algndrin,wlgndrin,meps,etrans2,successflag)
               if (.not.successflag) then
                 !write (58,*) i,j, ecut1, ecut2, meps, etrans2
                 algndrin = meps
                 wlgndrin = 0._dl
               end if !successflag
               do k=1,size(algndrin)
                 ein = algndrin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlgndrin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlgndrin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransii is large negative number
               ediff = etrans1 - etrans2
               !put in genin:
               !a and w frs \int_{\etransii}^{\etransi} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,etrans2,etrans1,successflag)
                 if (.not.successflag) then
                   !write (58,*) i,j, ecut1, ecut2, etrans1, etrans2
                   agenin = etrans2
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,etrans2,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\etransi} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransi is large negative number
               if (etrans1.gt.0._dl) then
               !a and w for \int_{\etransi}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans1,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = p1 - qout !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransi}^{\infty} dE_in
               end if !(etrans1.gt.0._dl)


             !\int_{\ecuti}^{\infty} dE_out:
             else if ((j.le.3*size(algndrout)+size(alagurout)).and.(j.ge.3*size(algndrout)+1)) then

               !if statement protects against machine imprecision,
               !set \int dE_in to GLag if \etransii is large negative number
               ediff = etrans2 - meps
               !put in genin:
               !a and w frs \int_{m_e}^{\etransii} dE_in
               if ((ediff.lt.200._dl).and.(ediff.ge.0._dl)) then !Use Legendre
                 agenin = in_bins_lgndr%abscissas
                 wgenin = in_bins_lgndr%weights
                 call cscgqf(1,agenin,wgenin,meps,etrans2,successflag)
                 if (.not.successflag) then
                   !write (58,*) i,j, ecut1, ecut2, meps, etrans2
                   agenin = meps 
                   wgenin = 0._dl
                 end if !successflag
               else
                 agenin = in_bins_lagur%abscissas
                 wgenin = in_bins_lagur%weights
                 call cscgqf(5,agenin,wgenin,meps,1._dl,successflag,0._dl)
                 wgenin = exp(in_bins_lagur%abscissas)*wgenin
               end if
               do k=1,size(agenin)
                 ein = agenin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = eout - p1 + ein + qin !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wgenin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wgenin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{m_e}^{\etransii} dE_in

               !if statement protects against machine imprecision,
               !set \int dE_in to zero if \etransii is large negative number
               if (etrans2.gt.0._dl) then
               !a and w for \int_{\etransii}^{\infty} dE_in
               alagurin = in_bins_lagur%abscissas
               wlagurin = in_bins_lagur%weights
               call cscgqf(5,alagurin,wlagurin,etrans2,1._dl,successflag,0._dl)
               wlagurin = exp(in_bins_lagur%abscissas)*wlagurin
               do k=1,size(alagurin)
                 ein = alagurin(k) !E_in energy
                 qin = sqrt(ein**2 - meps**2) !q_in momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,3) = fd_equil_calc(ein,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 p4 = eout + ein - p1 !p_4 energy/momentum (different energy cons. expression)
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 toplim = p1 + qout !top limit
                 botlim = eout - p1 + ein - qin !bottom limit
                 call calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,.true.) !value from \int dy
                 kinemfactin = wlagurin(k)*weightin
                 call get_f_nuepma(woccarray,.true.,net2,frs2)
                 call calc_weight_nuepma(weightinswap,toplim,botlim,p1,eout,meps,.false.) !value from \int dy
                 kinemfactinswap = wlagurin(k)*weightinswap
                 call get_f_nuepma(woccarray,.false.,net2swap,frs2swap)
                 net1 = net1 + net2*kinemfactin + net2swap*kinemfactinswap
                 frs1 = frs1 + frs2*kinemfactin + frs2swap*kinemfactinswap
                 !write (116,*) ein, weightin(1,1)
                 !write (119,*) ein, weightin(1,2)
                 !write (120,*) ein, weightin(1,1)*frs2(1,1) + weightinswap(1,1)*frs2swap(1,1)
                 !write (121,*) ein, weightin(2,1)*frs2(2,1) + weightinswap(2,1)*frs2swap(2,1)
                 !tempint = tempint + frs2(1,1)*kinemfactin(1,1) + frs2swap(1,1)*kinemfactinswap(1,1)
               end do !\int_{\etransii}^{\infty} dE_in
               end if !(etrans2.gt.0._dl)


             end if !\int dE_out


             !write (116,'(/)')
             !write (118,*) eout, tempint

             net(4,:,:,i) = net(4,:,:,i) + net1(:,:)*kinemfact1
             frs(4,:,:,i) = frs(4,:,:,i) + frs1(:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         end if


       end do

       !if (abs(meps-3.2732992268804195_dl).lt.0.01_dl) then
       !  call mpi_barrier(mpi_comm_world,ierr)
       !  if (rank.eq.0) write (*,*) 'here1'
         !write (200 + rank,*) rank, net(4,1,1,2), frs(4,1,1,2)
         !call flush(200 + rank)
         !if (rank.eq.0) write(*,*) out_bins_lagur%abscissas
       !end if

       net = prefac*net
       frs = prefac*frs

       !do j=1,100
       !  write (124,*) allarray(1,j,:)
       !  write (124,*) allarray(2,j,:)
       !  write (124,*) allarray(3,j,:)
       !  write (124,*)
       !end do



       !deallocate(algndrout)
       !deallocate(wlgndrout)
       !deallocate(alagurout)
       !deallocate(wlagurout)
       !deallocate(algndrin)
       !deallocate(wlgndrin)
       !deallocate(alagurin)
       !deallocate(wlagurin)
       !deallocate(agenout)
       !deallocate(wgenout)
       !deallocate(agenin)
       !deallocate(wgenin)
       !deallocate(eouta)
       !deallocate(eoutw)

       !do i=1,nbins
       !  !write (114,*) eps_bins%abscissas(i),((testsum(i,j) - 8._dl/3._dl)*3._dl/8._dl,j=1,nbins)
       !  write (117,*) eps_bins%abscissas(i),net(4,2,1,i),frs(4,2,1,i) &
       !                ,(net(4,2,1,i)/frs(4,2,1,i))
       !end do
       !write (117,*) (net(4,1,1,i)/frs(4,1,1,i),i=1,nbins)

       !deallocate(testsum)


       return


       end subroutine annihilation_nuepma


!-------------------------------------------------------------------     
       
       subroutine calc_weight_nuepma(weightin,toplim,botlim,p1,eout,meps,swapflag)


!------linkages.

!      called by - [subroutine] scattering_nuer2
!      calls     - [function] r2int_m_1, r2int_m_2


!------remarks.

!      Calculates weight function for R_2 integral for \nu-e^\pm scattering.


       implicit none


!------throughput variables.

       real(dl), dimension(:,:), intent(out) :: weightin !output array
       real(dl), intent(in) :: toplim !first input
       real(dl), intent(in) :: botlim !second input
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: eout !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity
       logical, intent(in) :: swapflag !Flag for whether E_out = E_2 (.true.) or E_out = E_3


!------local variables.

       integer m,n !indicies
       real(dl) flafact !factors
       real(dl) term1, term2 !terms


!------procedure.


       flafact = 1._dl
       do m=1,3 !neutrino flavor
         do n=1,2 !can't sum over E_2 and E_3.
           if ((swapflag.and.(n.eq.1)).or.(.not.swapflag.and.(n.eq.2))) then
             term1 = 32._dl*(2._dl*sin2thw + flafact)**2
             term2 = 2._dl*sin2thw/(2._dl*sin2thw + flafact)*meps**2 !change in overall sign
           else
             term1 = 128._dl*sin2thw**2
             term2 = (2._dl*sin2thw + flafact)/2._dl/sin2thw*meps**2
           end if
           weightin(m,n) = term1*( &
                                    rint_l_1(toplim,p1,eout,meps) - rint_l_1(botlim,p1,eout,meps) &
                           + term2*(rint_l_2(toplim,p1,eout,meps) - rint_l_2(botlim,p1,eout,meps)))
         end do !n
         flafact = -1._dl
       end do !m


       return


       end subroutine calc_weight_nuepma


!-------------------------------------------------------------------     
       
       real(dl) function rint_l_1(x,p1,eout,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuepma
!      calls     - none


!------remarks.

!      Calculates part of weight function for e^\pm annihilation.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: eout !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       rint_l_1 = 0.25_dl*(((p1-eout)**2 - meps**2)**2*x &
                         - 2._dl/3._dl*((p1-eout)**2 - meps**2)*x**3 &
                         + 0.2_dl*x**5)


       return


       end function rint_l_1


!-------------------------------------------------------------------     
       
       real(dl) function rint_l_2(x,p1,eout,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuepma
!      calls     - none


!------remarks.

!      Calculates part of weight function for e^\pm annihilation.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: eout !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       rint_l_2 = 0.5_dl*(1._dl/3._dl*x**3 - ((p1-eout)**2 - meps**2)*x)


       return


       end function rint_l_2


!-----------------------------------------------------------

       subroutine get_f_nuepma(woccarray,swapflag,net,frs)


!------linkages.

!      called by - [subroutine] annihilation_nuepma
!      calls     - none

!------remarks.

!      Subroutine to calculate F(p1,e2,e3,p4).
!      Identical to get_f_nuer1


       implicit none


!------throughput variables.

       real(dl), dimension(:,:,:), intent(in) :: woccarray !Array of occupation probs.
       logical, intent(in) :: swapflag !Flag for whether E_out = E_2 (.true.) or E_out = E_3
       real(dl), dimension(:,:), intent(out) :: net !net values
       real(dl), dimension(:,:), intent(out) :: frs !forward-reverse summed values


!------local variables.

       integer n,m,nn !indicies
       real(dl) forfact, revfact


!------procedure.


!------calculated values of frs and net rates.

       !24May2015 EG: Fixed bug with n and nn for particles 2 and 3
       do n=1,2
         nn = 2 - n/2
         do m=1,3 !\nu flavor
           if (swapflag) then
             forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(m,nn,4)) &
                       *woccarray(4,1,3)*woccarray(4,2,2)
             revfact = woccarray(m,n,1)*woccarray(m,nn,4) &
                       *(1._dl - woccarray(4,1,3))*(1._dl - woccarray(4,2,2))
           else !permute indicies 2 <-> 3
             forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(m,nn,4)) &
                       *woccarray(4,1,2)*woccarray(4,2,3)
             revfact = woccarray(m,n,1)*woccarray(m,nn,4) &
                       *(1._dl - woccarray(4,1,2))*(1._dl - woccarray(4,2,3))
           end if
           frs(m,n) = forfact + revfact
           net(m,n) = forfact - revfact
         end do !m
       end do !n


       return


       end subroutine get_f_nuepma


       end module trans_nuepma
