
       module trans_nuer1


!------linkages.

!      used by main

!------remarks.

!      Module for calculating R_1 integral of \nu + e -> e + \nu


       use mainvar
       use transvar
       use interp


       implicit none


       contains


!-----------------------------------------------------------

       !subroutine scattering_nuer1(nutrans,eps_bins,e2_bins_leg,e2_bins_lag,e3_bins_leg,net,frs,meps,tcmpl,phie)
       subroutine scattering_nuer1(nutrans,e2_bins_leg,e2_bins_lag,e3_bins_leg,net,frs,meps,tcmpl,phie)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] cscgqf, calc_weight_nuer1, interp_occfrac
!                               ,get_f_nuer1


!-----remarks.

!     Calculates scattering rates for R_1 integral of \nu scattering on e^\pm.


!------modules.

       use gauss


       implicit none


       save


!------throughput variables.

       type (nuvar), intent(in), dimension(:,:) :: nutrans !nu occupation fractions and info.
       !type (bin_scheme), intent(in) :: eps_bins !epsilon abscissas and weights
       type (bin_scheme), intent(in) :: e2_bins_leg !Gauss-Legendre abscissas and weights for e2 integral
       type (bin_scheme), intent(in) :: e2_bins_lag !Gauss-Laguerre abscissas and weights for e2 integral
       type (bin_scheme), intent(in) :: e3_bins_leg !Gauss-Legendre abscissas and weights for e3 integral
       real(dl), intent(out), dimension(:,:,:,:,:) :: net !net scattering rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:,:) :: frs !forward-reverse summed scattering rates/tnmev**5
       real(dl), intent(in) :: meps !m_e/T_{cm}
       real(dl), intent(in) :: tcmpl !T_{cm}/T_{pl}
       real(dl), intent(in) :: phie !electron degeneracy parameter


!------local variables.

       integer i,j,k,m,mm,n
       !integer nbins
       integer sumind

       real(dl) p1, e2, q2, e3, q3, p4
       real(dl) prefac
       real(dl) kinemfact1
       real(dl), dimension(2,3,2) :: kinemfact2
       real(dl) parfact
       real(dl) pcase1
       real(dl) ecut1, ecut3
       real(dl) etrans2
       real(dl) elim1, elim2

       real(dl), dimension(3,2) :: logwoccfrac

       real(dl), dimension(4,2,4) :: woccarray !array for calculating F

       logical successflag

       !real(dl), dimension(:), allocatable :: glegai2
       !real(dl), dimension(:), allocatable :: glegwi2

       !real(dl), dimension(:), allocatable :: glagai2
       !real(dl), dimension(:), allocatable :: glagwi2

       !real(dl), dimension(:), allocatable :: glegai3
       !real(dl), dimension(:), allocatable :: glegwi3

       !real(dl), dimension(:), allocatable :: e2a
       !real(dl), dimension(:), allocatable :: e2w

       !real(dl), dimension(:), allocatable :: e3a
       !real(dl), dimension(:), allocatable :: e3w

       real(dl), dimension(2,3,2) :: weight3

       !real(dl) forfact, revfact

       real(dl), dimension(2,3,2) :: net1
       real(dl), dimension(2,3,2) :: frs1
       real(dl), dimension(2,3,2) :: net2
       real(dl), dimension(2,3,2) :: frs2

       real(dl) botlim, toplim

       integer startind !starting index


!------Testing quantities.

       !real(dl) lowsum, highsum
       !real(dl), dimension(:,:), allocatable :: testsum
       !real(dl) tempint
       !real(dl) testsum2, woccfrac_exact


!------procedure.


       !allocate(glegai2(size(e2_bins_leg%abscissas))) !assign dimension of glegai2
       !allocate(glegwi2(size(e2_bins_leg%abscissas))) !assign dimension of glegwi2
       
       !allocate(glagai2(size(e2_bins_lag%abscissas))) !assign dimension of glagai2
       !allocate(glagwi2(size(e2_bins_lag%abscissas))) !assign dimension of glagwi2

       !allocate(glegai3(size(e3_bins_leg%abscissas))) !assign dimension of glegai3
       !allocate(glegwi3(size(e3_bins_leg%abscissas))) !assign dimension of glegwi3

       !allocate(e2a(2*size(glegai2)+size(glagai2))) !assign dimension of e2a
       !allocate(e2w(2*size(glegai2)+size(glagai2))) !assign dimension of e2w

       !allocate(e3a(size(glegai3))) !assign dimension of e3a
       !allocate(e3w(size(glegai3))) !assign dimension of e3w

       !nbins = size(nutrans(1,1)%occfrac)

       !allocate(testsum(nbins,size(e2_bins_leg%abscissas)+size(e2_bins_lag%abscissas))) !assign dimension of testsum

       do n=1,2
         do m=1,3 !don't use sterile occupation fractions
           woccfrac(m,n,:) = nutrans(m,n)%occfrac(:)
         end do
       end do
       woccarray = 0._dl


       !write (112,*) woccfrac(1,1,:)
       !do i=1,nbins
       !  write (113,*) eps_bins%abscissas(i), dlog(woccfrac(1,1,i)), woccfracddot(1,1,i)
       !end do

       prefac = gfhbar/16._dl/(2._dl*pi)**3
       !testsum = 0._dl
       !testsum2 = 0._dl
       !lowsum = 0._dl
       !highsum = 0._dl

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       pcase1 = meps/2._dl


       net = 0._dl
       frs = 0._dl
       sumind = 1
       do i=startind,nbins !outermost loop over p_1
       !do i=startind,nbins,10 !outermost loop over p_1
         p1 = eps_bins%abscissas(i)
         ecut1 = meps + 2._dl*p1**2/(meps - 2._dl*p1) !\ecuti, when \etransii is valid
         ecut3 = sqrt(p1**2 + meps**2) !\ecutiii, when q_2 = p_1
         woccarray(1:3,:,1) = woccfrac(1:3,:,i)

         !write (116,'(A,1pe12.5,/)') 'p_1/m_e = ', p1/meps
         !write (118,'(A,1pe12.5)') 'p1 = ', p1

         !p1/meps < 1/2:
         if (p1.lt.pcase1) then !ecut1 neccessary

           !Make the abscissas and weights for \int_{m_e}^{\infty} dE_2

           !put the input Gauss Legendre (GLeg) abscissas into a dummy array:
           glegai2 = e2_bins_leg%abscissas
           !put the input Gauss Legendre (GLeg) weights into a dummy array:
           glegwi2 = e2_bins_leg%weights
           !calculate the scaled abscissas and weights over given interval [meps,\ecutiii]
           call cscgqf(1,glegai2,glegwi2,meps,ecut3,successflag)
           !a and w for \int_{m_e}^{\ecutiii} dE_2
           e2a(1:size(glegai2)) = glegai2
           e2w(1:size(glegai2)) = glegwi2

           glegai2 = e2_bins_leg%abscissas
           glegwi2 = e2_bins_leg%weights
           call cscgqf(1,glegai2,glegwi2,ecut3,ecut1,successflag)
           !a and w for \int_{\ecutiii}^{\ecuti} dE_2
           e2a(size(glegai2)+1:2*size(glegai2)) = glegai2
           e2w(size(glegai2)+1:2*size(glegai2)) = glegwi2

           !Gauss Laguerre (GLag) interval:
           glagai2 = e2_bins_lag%abscissas
           glagwi2 = e2_bins_lag%weights
           call cscgqf(5,glagai2,glagwi2,ecut1,1._dl,successflag,0._dl)
           !a and w for \int_{\ecuti}^{\infty} dE_2
           e2a(2*size(glegai2)+1:2*size(glegai2)+size(glagai2)) = glagai2
           e2w(2*size(glegai2)+1:2*size(glegai2)+size(glagai2)) = exp(e2_bins_lag%abscissas)*glagwi2


           do j=1,2*size(glegai2)+size(glagai2)
           !do j=1,2*size(glegai2)+size(glagai2),15
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e2 = e2a(j)
             q2 = sqrt(e2**2 - meps**2)
             etrans2 = 0.5_dl*(2._dl*p1 + e2 - q2 + meps**2/(2._dl*p1 + e2 - q2))
             elim1 = 0.5_dl*(2._dl*p1 + e2 + q2 + meps**2/(2._dl*p1 + e2 + q2))
             elim2 = etrans2
             parfact = 1._dl
             do n=1,2
               !woccarray(4,n,2) = 1._dl/(exp(e2*tcmpl - parfact*phie) + 1._dl)
               woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e2w(j)/p1**2
             net1 = 0._dl
             frs1 = 0._dl
             !if (i.eq.50) write (122,*) ecut1, ecut3

             !write (116,'(A,1pe12.5)') 'e2 = ', e2
             !tempint = 0._dl

             !\int_{m_e}^{\ecutiii} dE_2:
             if (j.le.size(glegai2)) then

               !a and w for \int_{m_e}^{E_2} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,meps,e2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, meps, e2
                 glegai3 = meps
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k) !E_3 energy
                 q3 = sqrt(e3**2 - meps**2) !q_3 momentum
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3 !top limit
                 botlim = p1 + e2 - e3 - q3 !bottom limit
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps) !value from \int dy
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{m_e}^{E_2} dE_3

               !a and w for \int_{E_2}^{\etransii} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,e2,etrans2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, e2, etrans2
                 glegai3 = e2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = p1 - q2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{E_2}^{\etransii} dE_3

               !a and w for \int_{\etransii}^{\elimi} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,etrans2,elim1,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, etrans2, elim1
                 glegai3 = etrans2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = e3 + q3 - p1 - e2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{\etransii}^{\elimi} dE_3


             !\int_{\ecutiii}^{\ecuti} dE_2
             else if ((j.le.2*size(glegai2)).and.(j.ge.size(glegai2)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, meps, etrans2
                 glegai3 = meps
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k) !E_3 energy
                 q3 = sqrt(e3**2 - meps**2) !q_3 momentum
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3 !top limit
                 botlim = p1 + e2 - e3 - q3 !bottom limit
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 !woccarray(1:3,:,4) =  1._dl/(exp(p4) + 1._dl)
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (123,*) e3, net1(1,1,1), frs1(1,1,1), net1(1,1,1)/frs1(1,1,1), p4, woccarray(1,1,4)
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{m_e}^{\etransii} dE_3

               !a and w for \int_{\etransii}^{E_2} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,etrans2,e2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, etrans2, e2
                 glegai3 = etrans2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3
                 botlim = q2 - p1
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 !woccarray(1:3,:,4) =  1._dl/(exp(p4) + 1._dl)
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (123,*) e3, net1(1,1,1), frs1(1,1,1), net1(1,1,1)/frs1(1,1,1), p4, woccarray(1,1,4)
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{\etransii}^{E_2} dE_3

               !a and w for \int_{E_2}^{\elimi} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,e2,elim1,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, e2, elim1
                 glegai3 = e2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = e3 + q3 - p1 - e2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 !woccarray(1:3,:,4) =  1._dl/(exp(p4) + 1._dl)
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (123,*) e3, net1(1,1,1), frs1(1,1,1), net1(1,1,1)/frs1(1,1,1), p4, woccarray(1,1,4)
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{E_2}^{\elimi} dE_3


             !\int_{\ecuti}^{\infty} dE_2
             else

               !a and w for \int_{\elimii}^{E_2} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,elim2,e2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, elim2, e2
                 glegai3 = elim2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3
                 botlim = q2 - p1
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{\elimii}^{E_2}dE_3

               !a and w for \int_{E_2}^{\elimi} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,e2,elim1,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, ecut1, e2, elim1
                 glegai3 = e2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = e3 + q3 - p1 - e2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{E_2}^{\elimi} dE_3


             end if
             !write (123,'(/)')

             !write (116,'(/)')
             !write (118,*) e2, tempint

             !if (i.eq.50) then
             !  write (122,*) e2, e2w(j), net1(1,1,1), frs1(1,1,1), net1(1,1,1)/frs1(1,1,1)
             !end if
             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1
             !write (124,*) e2, e2w(j), net(4,1,1,1,50), frs(4,1,1,1,50), net(4,1,1,1,50)/frs(4,1,1,1,50)


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do !\int dE_2

           !write (118,'(/)')


         !1/2 < p1/meps:
         else !ecut1 unneccessary


           glegai2 = e2_bins_leg%abscissas
           glegwi2 = e2_bins_leg%weights
           call cscgqf(1,glegai2,glegwi2,meps,ecut3,successflag)
           !a and w for \int_{m_e}^{\ecutiii} dE_2
           e2a(1:size(glegai2)) = glegai2
           e2w(1:size(glegai2)) = glegwi2

           !Gauss Laguerre (GLag) interval:
           glagai2 = e2_bins_lag%abscissas
           glagwi2 = e2_bins_lag%weights
           call cscgqf(5,glagai2,glagwi2,ecut3,1._dl,successflag,0._dl)
           !a and w for \int_{\ecutiii}^{\infty} dE_2
           e2a(size(glegai2)+1:size(glegai2)+size(glagai2)) = glagai2
           e2w(size(glegai2)+1:size(glegai2)+size(glagai2)) = exp(e2_bins_lag%abscissas)*glagwi2


           do j=1,size(glegai2)+size(glagai2)
           !do j=1,size(glegai2)+size(glagai2),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e2 = e2a(j)
             q2 = sqrt(e2**2 - meps**2)
             etrans2 = 0.5_dl*(2._dl*p1 + e2 - q2 + meps**2/(2._dl*p1 + e2 - q2))
             elim1 = 0.5_dl*(2._dl*p1 + e2 + q2 + meps**2/(2._dl*p1 + e2 + q2))
             parfact = 1._dl
             do n=1,2
               !woccarray(4,n,2) = 1._dl/(exp(e2*tcmpl - parfact*phie) + 1._dl)
               woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e2w(j)/p1**2
             net1 = 0._dl
             frs1 = 0._dl
             !if (i.eq.50) write (122,*) p1 + e2, elim1

             !write (116,'(A,1pe12.5)') 'e2 = ', e2
             !tempint = 0._dl

             !\int_{m_e}^{\ecutiii} dE_2:
             if (j.le.size(glegai2)) then

               !a and w for \int_{m_e}^{E_2} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,meps,e2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, meps, e2
                 glegai3 = meps
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k) !E_3 energy
                 q3 = sqrt(e3**2 - meps**2) !q_3 momentum
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3 !top limit
                 botlim = p1 + e2 - e3 - q3 !bottom limit
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{m_e}^{E_2} dE_3

               !a and w for \int_{E_2}^{\etransii} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,e2,etrans2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, e2, etrans2
                 glegai3 = e2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = p1 - q2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{E_2}^{\etransii} dE_3

               !a and w for \int_{\etransii}^{\elimi} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,etrans2,elim1,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, etrans2, elim1
                 glegai3 = etrans2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = e3 + q3 - p1 - e2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{\etransii}^{\elimi} dE_3


             !\int_{\ecutiii}^{\infty} dE_2
             else

               !a and w for \int_{m_e}^{\etransii} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, meps, etrans2
                 glegai3 = meps
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k) !E_3 energy
                 q3 = sqrt(e3**2 - meps**2) !q_3 momentum
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3 !top limit
                 botlim = p1 + e2 - e3 - q3 !bottom limit
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{m_e}^{\etransii} dE_3

               !a and w for \int_{\etransii}^{E_2} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,etrans2,e2,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, etrans2, e2
                 glegai3 = etrans2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + e2 - e3 + q3
                 botlim = q2 - p1
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{\etransii}^{E_2} dE_3

               !a and w for \int_{E_2}^{\elimi} dE_3
               glegai3 = e3_bins_leg%abscissas
               glegwi3 = e3_bins_leg%weights
               call cscgqf(1,glegai3,glegwi3,e2,elim1,successflag)
               if (.not.successflag) then
                 write (59,*) i,j, ecut3, e2, elim1
                 glegai3 = e2
                 glegwi3 = 0._dl
               end if !successflag
               do k=1,size(glegai3)
                 e3 = glegai3(k)
                 q3 = sqrt(e3**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   !woccarray(4,n,3) = 1._dl/(exp(e3*tcmpl - parfact*phie) + 1._dl)
                   woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q2
                 botlim = e3 + q3 - p1 - e2
                 call calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)
                 kinemfact2 = glegwi3(k)*weight3
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer1(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !testsum(i,j) = testsum(i,j) + glegwi3(k)*weight3(1,1,1) &
                 !               /p1**3/e2**3/32._dl/(2._dl*sin2thw + 1._dl)**2
                 !write (116,*) e3, weight3(1,2,1)
                 !write (116,*) e3, weight3(1,2,1)*frs2(1,2,1)
                 !tempint = tempint + frs2(1,2,1)*kinemfact2(1,2,1)
               end do !\int_{E_2}^{\elimi} dE_3


             end if
             !if (i.eq.50) write (120,'(/)')

             !write (116,'(/)')
             !write (118,*) e2, tempint

             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do !\int dE_2

           !write (118,'(/)')


         end if !p_1 cases


       end do !p_1 loop

       net = prefac*net
       frs = prefac*frs


       !deallocate(glegai2)
       !deallocate(glegwi2)
       !deallocate(glagai2)
       !deallocate(glagwi2)
       !deallocate(glegai3)
       !deallocate(glegwi3)
       !deallocate(e2a)
       !deallocate(e2w)
       !deallocate(e3a)
       !deallocate(e3w)

       !do i=1,nbins
       !  !write (114,*) eps_bins%abscissas(i),((testsum(i,j) - 8._dl/3._dl)*3._dl/8._dl, &
       !  !                                     j=1,size(e2_bins_leg%abscissas)+size(e2_bins_lag%abscissas))
       !  write (117,*) eps_bins%abscissas(i),eps_bins%weights(i),net(4,1,1,1,i),frs(4,1,1,1,i) &
       !                ,(net(4,1,1,1,i)/frs(4,1,1,1,i))
       !end do
       !write (117,*) (net(4,1,1,1,i)/frs(4,1,1,1,i),i=1,nbins)


       !deallocate(testsum)


       return


       end subroutine scattering_nuer1


!-------------------------------------------------------------------     
       
       subroutine calc_weight_nuer1(weight3,toplim,botlim,p1,e2,meps)


!------linkages.

!      called by - [subroutine] scattering_nuer2
!      calls     - [function] r1int_m_1, r1int_m_2


!------remarks.

!      Calculates weight function for R_1 integral for \nu-e^\pm scattering.


       implicit none


!------throughput variables.

       real(dl), dimension(:,:,:), intent(out) :: weight3 !output array
       real(dl), intent(in) :: toplim !first input
       real(dl), intent(in) :: botlim !second input
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e2 !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------local variables.

       integer m,n,nn !indicies
       real(dl) flafact !factors
       real(dl) term1, term2 !terms


!------procedure.


       do n=1,2 !neutrino handedness
         flafact = 1._dl
         do m=1,3 !neutrino flavor
           do nn=1,2 !lepton type
             if (n.eq.nn) then
               term1 = 32._dl*(2._dl*sin2thw + flafact)**2
               term2 = -2._dl*sin2thw/(2._dl*sin2thw + flafact)*meps**2
             else
               term1 = 128._dl*sin2thw**2
               term2 = (2._dl*sin2thw + flafact)/2._dl/sin2thw*meps**2
             end if
             weight3(nn,m,n) = term1*( &
                                        r1int_m_1(toplim,p1,e2,meps) - r1int_m_1(botlim,p1,e2,meps) &
                               + term2*(r1int_m_2(toplim,p1,e2,meps) - r1int_m_2(botlim,p1,e2,meps)))
           end do !nn
           flafact = -1._dl
         end do !m
       end do !n


       return


       end subroutine calc_weight_nuer1


!-------------------------------------------------------------------     
       
       real(dl) function r1int_m_1(x,p1,e2,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuer1
!      calls     - none


!------remarks.

!      Calculates part of weight function for \nu-e scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e2 !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       r1int_m_1 = 0.25_dl*(((p1+e2)**2 - meps**2)**2*x &
                         - 2._dl/3._dl*((p1+e2)**2 - meps**2)*x**3 &
                         + 0.2_dl*x**5)


       return


       end function r1int_m_1


!-------------------------------------------------------------------     
       
       real(dl) function r1int_m_2(x,p1,e2,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuer1
!      calls     - none


!------remarks.

!      Calculates part of weight function for \nu-e scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e2 !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       r1int_m_2 = 0.5_dl*(((p1+e2)**2 - meps**2)*x - 1._dl/3._dl*x**3)


       return


       end function r1int_m_2


!-----------------------------------------------------------

       subroutine get_f_nuer1(woccarray,net,frs)


!------linkages.

!      called by - [program] scattering_nuer1
!      calls     - none

!------remarks.

!      Subroutine to calculate F(p1,e2,e3,p4).


       implicit none


!------throughput variables.

       real(dl), dimension(:,:,:), intent(in) :: woccarray !Array of occupation probs.
       real(dl), dimension(:,:,:), intent(out) :: net !net values
       real(dl), dimension(:,:,:), intent(out) :: frs !frs values


!------local variables.

       integer n,m,nn !indicies
       real(dl) forfact, revfact


!------procedure.


!------calculated values of frs and net rates.

       do n=1,2 !\nu handedness
         do m=1,3 !\nu flavor
           do nn=1,2 !lepton parity
             forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(4,nn,2)) &
                       *woccarray(4,nn,3)*woccarray(m,n,4)
             revfact = woccarray(m,n,1)*woccarray(4,nn,2) &
                       *(1._dl - woccarray(4,nn,3))*(1._dl - woccarray(m,n,4))
             frs(nn,m,n) = forfact + revfact
             net(nn,m,n) = forfact - revfact
           end do !nn
         end do !m
       end do !n


       return


       end subroutine get_f_nuer1


       end module trans_nuer1
