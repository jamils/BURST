
       module trans_nuer2


!------linkages.

!      used by main

!------remarks.

!      Module for calculating R_2 integral of \nu + e -> e + \nu


       use mainvar
       use transvar
       use interp


       implicit none


       contains


!-----------------------------------------------------------

       !subroutine scattering_nuer2(nutrans,eps_bins,e2_bins_leg,e2_bins_lag &
       !                            ,e3_bins_leg,e3_bins_lag,net,frs,meps,tcmpl,phie)
       subroutine scattering_nuer2(nutrans,e2_bins_leg,e2_bins_lag &
                                   ,e3_bins_leg,e3_bins_lag,net,frs,meps,tcmpl,phie)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] none
!                  [function] none


!-----remarks.

!     Calculates scattering rates for \nu scattering on e^\pm.


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
       type (bin_scheme), intent(in) :: e3_bins_lag !Gauss-Laguerre abscissas and weights for e3 integral
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
       real(dl) pcase1, pcase2, pcase3
       real(dl) ecut1, ecut2, ecut3
       real(dl) etrans2
       real(dl) elim1, elim2

       real(dl), dimension(3,2) :: logwoccfrac

       real(dl), dimension(4,2,4) :: woccarray !array for calculating F

       logical successflag, ecutflag

       !real(dl), dimension(:), allocatable :: glegai2
       !real(dl), dimension(:), allocatable :: glegwi2

       !real(dl), dimension(:), allocatable :: glagai2
       !real(dl), dimension(:), allocatable :: glagwi2

       !real(dl), dimension(:), allocatable :: glegai3
       !real(dl), dimension(:), allocatable :: glegwi3

       !real(dl), dimension(:), allocatable :: glagai3
       !real(dl), dimension(:), allocatable :: glagwi3

       !real(dl), dimension(:), allocatable :: ggenai2
       !real(dl), dimension(:), allocatable :: ggenwi2

       !real(dl), dimension(:), allocatable :: e3r2a
       !real(dl), dimension(:), allocatable :: e3r2w

       real(dl), dimension(2,3,2) :: weight2

       real(dl), dimension(2,3,2) :: net1
       real(dl), dimension(2,3,2) :: frs1
       real(dl), dimension(2,3,2) :: net2
       real(dl), dimension(2,3,2) :: frs2

       real(dl) botlim, toplim

       integer startind !starting index


!------Testing quantities.

       !real(dl), dimension(:,:), allocatable :: testsum
       !real(dl) tempint
       !real(dl) lowsum, highsum
       !real(dl) testsum2, woccfrac_exact
       !real(dl), dimension(3,100,3) :: allarray
       !integer jj


!------procedure.


       !allocate(glegai2(size(e2_bins_leg%abscissas))) !assign dimension of glegai2
       !allocate(glegwi2(size(e2_bins_leg%abscissas))) !assign dimension of glegwi2
       
       !allocate(glagai2(size(e2_bins_lag%abscissas))) !assign dimension of glagai2
       !allocate(glagwi2(size(e2_bins_lag%abscissas))) !assign dimension of glagwi2

       !allocate(glegai3(size(e3_bins_leg%abscissas))) !assign dimension of glegai3
       !allocate(glegwi3(size(e3_bins_leg%abscissas))) !assign dimension of glegwi3

       !allocate(glagai3(size(e3_bins_lag%abscissas))) !assign dimension of glagai3
       !allocate(glagwi3(size(e3_bins_lag%abscissas))) !assign dimension of glagwi3

       !allocate(ggenai2(max(size(e2_bins_lag%abscissas),size(e2_bins_leg%abscissas)))) !assign dimension of ggenai2
       !allocate(ggenwi2(max(size(e2_bins_lag%abscissas),size(e2_bins_leg%abscissas)))) !assign dimension of ggenwi2

       !allocate(e3r2a(3*size(glegai3)+size(glagai3))) !assign dimension of e3r2a
       !allocate(e3r2w(3*size(glegai3)+size(glagai3))) !assign dimension of e3r2w


       !nbins = size(nutrans(1,1)%occfrac)

       !allocate(testsum(nbins,3*size(e3_bins_leg%abscissas)+size(e3_bins_lag%abscissas))) !assign dimension of testsum
       


       do n=1,2
         do m=1,3 !don't use sterile occupation fractions
           woccfrac(m,n,:) = nutrans(m,n)%occfrac(:)
         end do
       end do
       woccarray = 0._dl


       prefac = gfhbar/16._dl/(2._dl*pi)**3
       !testsum = 0._dl
       !testsum2 = 0._dl

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       pcase1 = meps*(sqrt(5._dl) - 1._dl)/4._dl
       pcase2 = meps/2._dl/sqrt(2._dl)
       pcase3 = meps/2._dl


       net = 0._dl
       frs = 0._dl
       sumind = 1
       !write (116,'(A,1pe12.5)') 'm_e = ', meps
       !write (118,'(A,1pe12.5)') 'm_e = ', meps
       do i=startind,nbins !outermost loop over p_1
       !do i=startind,nbins,10 !outermost loop over p_1
       !do i=3,3 !outermost loop over p_1
         p1 = eps_bins%abscissas(i)
         ecut1 = p1 + meps**2/4._dl/p1 !=\ecuti, when \emaxb finite
         ecut2 = p1 + meps*(p1 + meps)/(2._dl*p1 + meps) !=\ecutii, when \etransii defined
         ecut3 = sqrt(p1**2 + meps**2) !=\ecutiii, when q_3 = p_1
         woccarray(1:3,:,1) = woccfrac(1:3,:,i)

         !write (116,*)
         !write (116,'(A,1pe12.5)') 'p_1 = ', p1
         !write (116,*) ecut3, ecut2, ecut1
         !write (118,*)
         !write (118,'(A,1pe12.5)') 'p1 = ', p1

         !p1/meps < (\sqrt{5}-1)/4:
         if (p1.lt.pcase1) then !ecut3 < ecut2 < ecut1

           !Make the abscissas and weights for \int_{m_e}^{\infty} dE_3

           !put the input Gauss Legendre (GLeg) abscissas into a dummy array:
           glegai3 = e3_bins_leg%abscissas
           !put the input Gauss Legendre (GLeg) weights into a dummy array:
           glegwi3 = e3_bins_leg%weights
           !calculate the scaled abscissas and weights over given interval [meps,\ecutiii]
           call cscgqf(1,glegai3,glegwi3,meps,ecut3,successflag)
           !a and w for \int_{m_e}^{\ecutiii} dE_3
           e3r2a(1:size(glegai3)) = glegai3
           e3r2w(1:size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut3,ecut2,successflag)
           !a and w for \int_{\ecutiii}^{\ecutii} dE_3
           e3r2a(size(glegai3)+1:2*size(glegai3)) = glegai3
           e3r2w(size(glegai3)+1:2*size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut2,ecut1,successflag)
           !a and w for \int_{\ecutii}^{\ecuti} dE_3
           e3r2a(2*size(glegai3)+1:3*size(glegai3)) = glegai3
           e3r2w(2*size(glegai3)+1:3*size(glegai3)) = glegwi3

           !Gauss Laguerre (GLag) interval:
           glagai3 = e3_bins_lag%abscissas
           glagwi3 = e3_bins_lag%weights
           call cscgqf(5,glagai3,glagwi3,ecut1,1._dl,successflag,0._dl)
           !a and w for \int_{\ecuti}^{\infty} dE_3
           e3r2a(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = glagai3
           e3r2w(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = exp(e3_bins_lag%abscissas)*glagwi3

           do k=1,3*size(glegai3)+size(glagai3)
           !do k=10,3*size(glegai3)+size(glagai3),10
           !jj = 0
           !do k=149,150
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e3 = e3r2a(k)
             q3 = sqrt(e3**2 - meps**2)
             etrans2 = 0.5_dl*(e3 + q3 - 2._dl*p1 + meps**2/(e3 + q3 - 2._dl*p1))
             !18Jan2016 EG: A slight issue that elim1 could be negative
             !due to machine precision:
             elim1 = 0.5_dl*(e3 - q3 - 2._dl*p1 + meps**2/(e3 - q3 - 2._dl*p1))
             if (elim1.lt.0._dl) elim1 = max(e3, etrans2)
             elim2 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e3r2w(k)/p1**2
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'e3 = ', e3
             !write (116,*) 'e3 = ', e3, etrans2, elim1, k
             !tempint = 0._dl
             !jj = jj + 1
             !write (*,*) elim2, e3, elim1

             !\int_{m_e}^{\ecutiii} dE_3:
             if (k.le.size(glegai3)) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{E_3} dE_2

               !a and w for \int_{E_3}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,e3,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, e3, etrans2
                 glegai2 = e3
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{\elimi} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,elim1,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, etrans2, elim1
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{\elimi} dE_2


             !\int_{\ecutiii}^{\ecutii} dE_3
             else if ((k.le.2*size(glegai3)).and.(k.ge.size(glegai3)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, meps, etrans2
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, etrans2, e3
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{E_3} dE_2

               !a and w for \int_{E_3}^{\elimi} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,e3,elim1,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, e3, elim1
                 glegai2 = e3
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\elimi} dE_2


             !\int_{\ecutii}^{\ecuti} dE_3
             else if ((k.le.3*size(glegai3)).and.(k.ge.2*size(glegai3)+1)) then

               !a and w for \int_{\elimii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,elim2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut2, ecut1, elim2, e3
                 glegai2 = elim2
                 glegwi2 = 0._dl
               end if
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(jj,j,1) = glegai2(j)
                 !allarray(jj,j,2) = glegwi2(j)
                 !allarray(jj,j,3) = frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\elimii}^{E_3} dE_2

               !a and w for \int_{E_3}^{\elimi} dE_2
               !Depending on E_3 and \elimi, use different GQ methods
               if ((elim1-e3).lt.200._dl) then !Use Legendre
                 ggenai2 = e2_bins_leg%abscissas
                 ggenwi2 = e2_bins_leg%weights
                 call cscgqf(1,ggenai2,ggenwi2,e3,elim1,successflag)
                 if (.not.successflag) then
                   write (60,*) i,k, ecut2, ecut1, e3, elim1
                   ggenai2 = e3
                   ggenwi2 = 0._dl
                 end if
               else !Use Laguerre
                 ggenai2 = e2_bins_lag%abscissas
                 ggenwi2 = e2_bins_lag%weights
                 call cscgqf(5,ggenai2,ggenwi2,e3,1._dl,successflag,0._dl)
                 ggenwi2 = exp(e2_bins_lag%abscissas)*ggenwi2
               end if
               do j=1,size(ggenai2)
                 e2 = ggenai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = ggenwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(jj,j+50,1) = ggenai2(j)
                 !allarray(jj,j+50,2) = ggenwi2(j)
                 !allarray(jj,j+50,3) = frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\elimi} dE_2


             !\int_{\ecuti}^{\infty} dE_3
             else

               !a and w for \int_{\elimii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,elim2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, elim2, e3
                 glegai2 = elim2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\elimii}^{E_3}dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             end if !\int dE_3


             !write (116,'(/)')
             !write (118,*) e3, tempint

             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !(\sqrt{5}-1)/4 < p1/meps < 1/(2*\sqrt{2}):
         else if ((p1.lt.pcase2).and.(p1.ge.pcase1)) then !ecut3 < ecut1 < ecut2

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,meps,ecut3,successflag)
           !a and w for \int_{m_e}^{\ecutiii} dE_3
           e3r2a(1:size(glegai3)) = glegai3
           e3r2w(1:size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut3,ecut1,successflag)
           !a and w for \int_{\ecutiii}^{\ecuti} dE_3
           e3r2a(size(glegai3)+1:2*size(glegai3)) = glegai3
           e3r2w(size(glegai3)+1:2*size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut1,ecut2,successflag)
           !a and w for \int_{\ecuti}^{\ecutii} dE_3
           e3r2a(2*size(glegai3)+1:3*size(glegai3)) = glegai3
           e3r2w(2*size(glegai3)+1:3*size(glegai3)) = glegwi3

           glagai3 = e3_bins_lag%abscissas
           glagwi3 = e3_bins_lag%weights
           call cscgqf(5,glagai3,glagwi3,ecut2,1._dl,successflag,0._dl)
           !a and w for \int_{\ecutii}^{\infty} dE_3
           e3r2a(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = glagai3
           e3r2w(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = exp(e3_bins_lag%abscissas)*glagwi3

           do k=1,3*size(glegai3)+size(glagai3)
           !do k=10,3*size(glegai3)+size(glagai3),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e3 = e3r2a(k)
             q3 = sqrt(e3**2 - meps**2)
             etrans2 = 0.5_dl*(e3 + q3 - 2._dl*p1 + meps**2/(e3 + q3 - 2._dl*p1))
             !18Jan2016 EG: A slight issue that elim1 could be negative
             !due to machine precision:
             elim1 = 0.5_dl*(e3 - q3 - 2._dl*p1 + meps**2/(e3 - q3 - 2._dl*p1))
             if (elim1.lt.0._dl) elim1 = max(e3, etrans2)
             elim2 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e3r2w(k)/p1**2
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'e3 = ', e3
             !tempint = 0._dl

             !\int_{m_e}^{\ecutiii} dE_3:
             if (k.le.size(glegai3)) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{E_3} dE_2

               !a and w for \int_{E_3}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,e3,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, e3, etrans2
                 glegai2 = e3
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{\elimi} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,elim1,successflag)
               if (.not.successflag) then
                 write (60,*) i,k,meps, ecut3, etrans2, elim1
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{\elimi} dE_2


             !\int_{\ecutiii}^{\ecuti} dE_3
             else if ((k.le.2*size(glegai3)).and.(k.ge.size(glegai3)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut1, meps, etrans2
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut1, etrans2, e3
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{E_3} dE_2

               !put in ggen:
               !a and w for \int_{E_3}^{\elimi} dE_2
               if ((elim1-e3).lt.200._dl) then !Use Legendre
                 ggenai2 = e2_bins_leg%abscissas
                 ggenwi2 = e2_bins_leg%weights
                 call cscgqf(1,ggenai2,ggenwi2,e3,elim1,successflag)
                 if (.not.successflag) then
                   write (60,*) i,k, ecut3, ecut1, e3, elim1
                   ggenai2 = e3
                   ggenwi2 = 0._dl
                 end if !successflag
               else !Use Laguerre
                 ggenai2 = e2_bins_lag%abscissas
                 ggenwi2 = e2_bins_lag%weights
                 call cscgqf(5,ggenai2,ggenwi2,e3,1._dl,successflag,0._dl)
                 ggenwi2 = exp(e2_bins_lag%abscissas)*ggenwi2
               end if
               do j=1,size(ggenai2)
                 e2 = ggenai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = ggenwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\elimi} dE_2


             !\int_{\ecuti}^{\ecutii} dE_3
             else if ((k.le.3*size(glegai3)).and.(k.ge.2*size(glegai3)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, ecut2, meps, etrans2
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, ecut2, etrans2, e3
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{E_3} dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             !\int_{\ecutii}^{\infty} dE_3
             else

               !a and w for \int_{\elimii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,elim2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut2, elim2, e3
                 glegai2 = elim2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\elimii}^{E_3}} dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             end if !\int dE_3

             !write (116,'(/)')
             !write (118,*) e3, tempint

             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !1/(2\sqrt{2}) < p1/meps < 1/2:
         else if ((p1.lt.pcase3).and.(p1.ge.pcase2)) then !ecut1 < ecut3 < ecut2
                                                          !etrans2 defined for E_3 < ecut1


           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,meps,ecut1,successflag)
           !a and w for \int_{m_e}^{\ecuti} dE_3
           e3r2a(1:size(glegai3)) = glegai3
           e3r2w(1:size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut1,ecut3,successflag)
           !a and w for \int_{\ecuti}^{\ecutiii} dE_3
           e3r2a(size(glegai3)+1:2*size(glegai3)) = glegai3
           e3r2w(size(glegai3)+1:2*size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut3,ecut2,successflag)
           !a and w for \int_{\ecutiii}^{\ecutii} dE_3
           e3r2a(2*size(glegai3)+1:3*size(glegai3)) = glegai3
           e3r2w(2*size(glegai3)+1:3*size(glegai3)) = glegwi3

           glagai3 = e3_bins_lag%abscissas
           glagwi3 = e3_bins_lag%weights
           call cscgqf(5,glagai3,glagwi3,ecut2,1._dl,successflag,0._dl)
           !a and w for \int_{\ecutii}^{\infty} dE_3
           e3r2a(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = glagai3
           e3r2w(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = exp(e3_bins_lag%abscissas)*glagwi3

           do k=1,3*size(glegai3)+size(glagai3)
           !do k=10,3*size(glegai3)+size(glagai3),10
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e3 = e3r2a(k)
             q3 = sqrt(e3**2 - meps**2)
             etrans2 = 0.5_dl*(e3 + q3 - 2._dl*p1 + meps**2/(e3 + q3 - 2._dl*p1))
             !18Jan2016 EG: A slight issue that elim1 could be negative
             !due to machine precision:
             elim1 = 0.5_dl*(e3 - q3 - 2._dl*p1 + meps**2/(e3 - q3 - 2._dl*p1))
             if (elim1.lt.0._dl) elim1 = max(e3, etrans2)
             elim2 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e3r2w(k)/p1**2
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'e3 = ', e3
             !tempint = 0._dl

             !\int_{m_e}^{\ecuti} dE_3:
             if (k.le.size(glegai3)) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, meps, ecut1, meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{E_3} dE_2

               !a and w for \int_{E_3}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,e3,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, meps, ecut1, e3, etrans2
                 glegai2 = e3
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{\elimi} dE_2
               if ((elim1-etrans2).lt.200._dl) then !Use Legendre
                 ggenai2 = e2_bins_leg%abscissas
                 ggenwi2 = e2_bins_leg%weights
                 call cscgqf(1,ggenai2,ggenwi2,etrans2,elim1,successflag)
                 if (.not.successflag) then
                   write (60,*) i,k, meps, ecut1, etrans2, elim1
                   ggenai2 = etrans2
                   ggenwi2 = 0._dl
                 end if !successflag
               else !Use Laguerre
                 ggenai2 = e2_bins_lag%abscissas
                 ggenwi2 = e2_bins_lag%weights
                 call cscgqf(5,ggenai2,ggenwi2,etrans2,1._dl,successflag,0._dl)
                 ggenwi2 = exp(e2_bins_lag%abscissas)*ggenwi2
               end if
               do j=1,size(ggenai2)
                 e2 = ggenai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = ggenwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{\elimi} dE_2


             !\int_{\ecuti}^{\ecutiii} dE_3
             else if ((k.le.2*size(glegai3)).and.(k.ge.size(glegai3)+1)) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, ecut3, meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransiii} dE_2

               !a and w for \int_{E_3}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,e3,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, ecut3, e3, etrans2
                 glegai2 = e3
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,etrans2,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{\infty} dE_2


             !\int_{\ecutiii}^{\ecutii} dE_3
             else if ((k.le.3*size(glegai3)).and.(k.ge.2*size(glegai3)+1)) then

               !a and w for \int_{m_e}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, meps, etrans2
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, etrans2, e3
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{E_3} dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             !\int_{\ecutii}^{\infty} dE_3
             else

               !a and w for \int_{\elimii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,elim2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut2, elim2, e3
                 glegai2 = elim2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\elimii}^{E_3} dE_2


               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             end if !\int dE_3

             !write (116,'(/)')
             !write (118,*) e3, tempint

             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         !1/2 < p1/meps:
         else !ecut1 < ecut3 < ecut2 and etrans2 not defined for E_3 < ecut1

           ecutflag = .true.

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,meps,ecut1,successflag)
           !a and w for \int_{m_e}^{\ecuti} dE_3
           e3r2a(1:size(glegai3)) = glegai3
           e3r2w(1:size(glegai3)) = glegwi3

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut1,ecut3,successflag)
           !a and w for \int_{\ecuti}^{\ecutiii} dE_3
           if (successflag) then
             e3r2a(size(glegai3)+1:2*size(glegai3)) = glegai3
             e3r2w(size(glegai3)+1:2*size(glegai3)) = glegwi3
           else !30Mar2015 EG: included for completeness in case m_e = 0; never used
             e3r2a(size(glegai3)+1:2*size(glegai3)) = ecut1
             e3r2w(size(glegai3)+1:2*size(glegai3)) = 0._dl
             ecutflag = .false.
           end if

           glegai3 = e3_bins_leg%abscissas
           glegwi3 = e3_bins_leg%weights
           call cscgqf(1,glegai3,glegwi3,ecut3,ecut2,successflag)
           !a and w for \int_{\ecutiii}^{\ecutii} dE_3
           if (successflag) then
             e3r2a(2*size(glegai3)+1:3*size(glegai3)) = glegai3
             e3r2w(2*size(glegai3)+1:3*size(glegai3)) = glegwi3
           else
             e3r2a(2*size(glegai3)+1:3*size(glegai3)) = ecut3
             e3r2w(2*size(glegai3)+1:3*size(glegai3)) = 0._dl
             ecutflag = .false.
           end if

           glagai3 = e3_bins_lag%abscissas
           glagwi3 = e3_bins_lag%weights
           call cscgqf(5,glagai3,glagwi3,ecut2,1._dl,successflag,0._dl)
           !a and w for \int_{\ecutii}^{\infty} dE_3
           e3r2a(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = glagai3
           e3r2w(3*size(glegai3)+1:3*size(glegai3)+size(glagai3)) = exp(e3_bins_lag%abscissas)*glagwi3

           !jj = 0
           do k=1,3*size(glegai3)+size(glagai3)
           !do k=10,3*size(glegai3)+size(glagai3),10
           !do k=49,51
#ifdef prllel
           if (rank.eq.mod(sumind,nprocs)) then
#endif
             e3 = e3r2a(k)
             q3 = sqrt(e3**2 - meps**2)
             if (meps.ne.0._dl) then
               etrans2 = 0.5_dl*(e3 + q3 - 2._dl*p1 + meps**2/(e3 + q3 - 2._dl*p1))
             else 
               etrans2 = e3 - p1
             end if
             elim2 = etrans2
             parfact = 1._dl
             do n=1,2
               woccarray(4,n,3) = fd_equil_calc(e3,tcmpl,parfact*phie)
               parfact = -1._dl
             end do
             kinemfact1 = e3r2w(k)/p1**2
             net1 = 0._dl
             frs1 = 0._dl

             !write (116,*)
             !write (116,*) 'e3 = ', e3
             !tempint = 0._dl
             !jj = jj + 1

             !\int_{m_e}^{\ecuti} dE_3:
             if (k.le.size(glegai3)) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, meps, ecut1, meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(jj,j,1) = frs2(1,1,1)
                 !allarray(jj,j,2) = kinemfact2(1,1,1)
                 !allarray(jj,j,3) = tempint
               end do !\int_{m_e}^{E_3} dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(jj,j+50,1) = frs2(1,1,1)
                 !allarray(jj,j+50,2) = kinemfact2(1,1,1)
                 !allarray(jj,j+50,3) = tempint
               end do !\int_{E_3}^{\infty} dE_2

               !write (*,*) e3, tempint


             !\int_{\ecuti}^{\ecutiii} dE_3
             else if ((k.le.2*size(glegai3)).and.(k.ge.size(glegai3)+1) &
                       .and.ecutflag) then

               !a and w for \int_{m_e}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut1, ecut3,  meps, e3
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransiii} dE_2

               !a and w for \int_{E_3}^{\etransii} dE_2
               if ((etrans2-e3).lt.200._dl) then !Use Legendre
                 ggenai2 = e2_bins_leg%abscissas
                 ggenwi2 = e2_bins_leg%weights
                 call cscgqf(1,ggenai2,ggenwi2,e3,etrans2,successflag)
                 if (.not.successflag) then
                   write (60,*) i,k, ecut1, ecut3, e3, etrans2
                   ggenai2 = e3
                   ggenwi2 = 0._dl
                 end if !successflag
               else !Use Laguerre
                 ggenai2 = e2_bins_lag%abscissas
                 ggenwi2 = e2_bins_lag%weights
                 call cscgqf(5,ggenai2,ggenwi2,e3,1._dl,successflag,0._dl)
                 ggenwi2 = exp(e2_bins_lag%abscissas)*ggenwi2
               end if
               do j=1,size(ggenai2)
                 e2 = ggenai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = p1 - q3
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = ggenwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,etrans2,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{\infty} dE_2


             !\int_{\ecutiii}^{\ecutii} dE_3
             else if ((k.le.3*size(glegai3)).and.(k.ge.2*size(glegai3)+1) &
                       .and.ecutflag) then

               !a and w for \int_{m_e}^{\etransii} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,meps,etrans2,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, meps, etrans2
                 glegai2 = meps
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j) !E_2 energy
                 q2 = sqrt(e2**2 - meps**2) !q_2 momentum
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2 !top limit
                 botlim = p1 - e3 + e2 - q2 !bottom limit
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3 !p_4 energy/momentum
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{m_e}^{\etransii} dE_2

               !a and w for \int_{\etransii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,etrans2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut3, ecut2, etrans2, e3
                 glegai2 = etrans2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{\etransii}^{E_3} dE_2

               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
               end do !\int_{E_3}^{\infty} dE_2


             !\int_{\ecutii}^{\infty} dE_3
             else if ((k.le.3*size(glegai3)+size(glagai3)).and.(k.ge.3*size(glegai3)+1)) then

               !a and w for \int_{\elimii}^{E_3} dE_2
               glegai2 = e2_bins_leg%abscissas
               glegwi2 = e2_bins_leg%weights
               call cscgqf(1,glegai2,glegwi2,elim2,e3,successflag)
               if (.not.successflag) then
                 write (60,*) i,k, ecut2, elim2, e3
                 glegai2 = elim2
                 glegwi2 = 0._dl
               end if !successflag
               do j=1,size(glegai2)
                 e2 = glegai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 - e3 + e2 + q2
                 botlim = q3 - p1
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glegwi2(j)*weight2
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(2,j,1) = frs2(1,1,1)
                 !allarray(2,j,2) = kinemfact2(1,1,1)
                 !allarray(2,j,3) = tempint
               end do !\int_{\elimii}^{E_3} dE_2


               !a and w for \int_{E_3}^{\infty} dE_2
               glagai2 = e2_bins_lag%abscissas
               glagwi2 = e2_bins_lag%weights
               call cscgqf(5,glagai2,glagwi2,e3,1._dl,successflag,0._dl)
               glagwi2 = exp(e2_bins_lag%abscissas)*glagwi2
               do j=1,size(glagai2)
                 e2 = glagai2(j)
                 q2 = sqrt(e2**2 - meps**2)
                 parfact = 1._dl
                 do n=1,2
                   woccarray(4,n,2) = fd_equil_calc(e2,tcmpl,parfact*phie)
                   parfact = -1._dl
                 end do
                 toplim = p1 + q3
                 botlim = e3 - p1 - e2 + q2
                 call calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps) !value from \int dy
                 kinemfact2 = glagwi2(j)*weight2 !01Apr2015 EG: there was a bug here; the weight was Legendre
                 p4 = p1 + e2 - e3
                 call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
                 woccarray(1:3,:,4) = dexp(logwoccfrac(1:3,:))
                 call get_f_nuer2(woccarray,net2,frs2)
                 net1 = net1 + net2*kinemfact2
                 frs1 = frs1 + frs2*kinemfact2
                 !write (116,*) e2, weight2(1,1,1)
                 !write (116,*) e2, weight2(1,1,1)*frs2(1,1,1)
                 !tempint = tempint + frs2(1,1,1)*kinemfact2(1,1,1)
                 !allarray(2,j+50,1) = frs2(1,1,1)
                 !allarray(2,j+50,2) = kinemfact2(1,1,1)
                 !allarray(2,j+50,3) = tempint
               end do !\int_{E_3}^{\infty} dE_2

               !write (*,*) e3, tempint


             end if !\int dE_3

             !write (116,'(/)')
             !if ((k.le.size(glegai2)).or.(k.gt.3*size(glegai2))) then
             !  write (118,*) e3, tempint
             !else
             !  if (ecutflag) write (118,*) e3, tempint
             !end if

             net(4,:,:,:,i) = net(4,:,:,:,i) + net1(:,:,:)*kinemfact1
             frs(4,:,:,:,i) = frs(4,:,:,:,i) + frs1(:,:,:)*kinemfact1


#ifdef prllel
           end if
#endif
             sumind = sumind + 1

           end do


         end if


       end do

       net = prefac*net
       frs = prefac*frs

       !do j=1,100
       !  write (124,*) allarray(1,j,:)
       !  write (124,*) allarray(2,j,:)
       !  write (124,*) allarray(3,j,:)
       !  write (124,*)
       !end do


       !deallocate(glegai2)
       !deallocate(glegwi2)
       !deallocate(glagai2)
       !deallocate(glagwi2)
       !deallocate(glegai3)
       !deallocate(glegwi3)
       !deallocate(glagai3)
       !deallocate(glagwi3)
       !deallocate(ggenai2)
       !deallocate(ggenwi2)
       !deallocate(e3r2a)
       !deallocate(e3r2w)

       !do i=1,nbins
         !write (114,*) eps_bins%abscissas(i),((testsum(i,j) - 8._dl/3._dl)*3._dl/8._dl,j=1,nbins)
         !write (117,*) eps_bins%abscissas(i),eps_bins%weights(i),net(1,1,1,1,i),frs(1,1,1,1,i) &
         !              ,(net(1,1,1,1,i)/frs(1,1,1,1,i))
       !end do
       !write (117,*) (net(1,1,1,1,i)/frs(1,1,1,1,i),i=1,nbins)

       !deallocate(testsum)

       !if (rank.eq.42) write (*,*) tcmpl, meps, net(4,:,1,1,9)

       return


       end subroutine scattering_nuer2


!-------------------------------------------------------------------     
       
       subroutine calc_weight_nuer2(weight2,toplim,botlim,p1,e3,meps)


!------linkages.

!      called by - [subroutine] scattering_nuer2
!      calls     - [function] r2int_m_1, r2int_m_2


!------remarks.

!      Calculates weight function for R_2 integral for \nu-e^\pm scattering.


       implicit none


!------throughput variables.

       real(dl), dimension(:,:,:), intent(out) :: weight2 !output array
       real(dl), intent(in) :: toplim !first input
       real(dl), intent(in) :: botlim !second input
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e3 !second epsilon quantity
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
               term1 = 128._dl*sin2thw**2
               term2 = (2._dl*sin2thw + flafact)/2._dl/sin2thw*meps**2
             else
               term1 = 32._dl*(2._dl*sin2thw + flafact)**2
               term2 = -2._dl*sin2thw/(2._dl*sin2thw + flafact)*meps**2
             end if
             weight2(nn,m,n) = term1*( &
                                        r2int_m_1(toplim,p1,e3,meps) - r2int_m_1(botlim,p1,e3,meps) &
                               + term2*(r2int_m_2(toplim,p1,e3,meps) - r2int_m_2(botlim,p1,e3,meps)))
           end do !nn
           flafact = -1._dl
         end do !m
       end do !n


       return


       end subroutine calc_weight_nuer2


!-------------------------------------------------------------------     
       
       real(dl) function r2int_m_1(x,p1,e3,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuer1
!      calls     - none


!------remarks.

!      Calculates part of weight function for \nu-e scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e3 !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       r2int_m_1 = 0.25_dl*(((p1-e3)**2 - meps**2)**2*x &
                         - 2._dl/3._dl*((p1-e3)**2 - meps**2)*x**3 &
                         + 0.2_dl*x**5)


       return


       end function r2int_m_1


!-------------------------------------------------------------------     
       
       real(dl) function r2int_m_2(x,p1,e3,meps)


!------linkages.

!      called by - [subroutine] calc_weight_nuer1
!      calls     - none


!------remarks.

!      Calculates part of weight function for \nu-e scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !input limit quantity
       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: e3 !second epsilon quantity
       real(dl), intent(in) :: meps !m_e epsilon quantity


!------procedure.


       r2int_m_2 = 0.5_dl*(1._dl/3._dl*x**3 - ((p1-e3)**2 - meps**2)*x)


       return


       end function r2int_m_2


!-----------------------------------------------------------

       subroutine get_f_nuer2(woccarray,net,frs)


!------linkages.

!      called by - [program] scattering_nuer2
!      calls     - none

!------remarks.

!      Subroutine to calculate F(p1,e2,e3,p4).
!      Identical to get_f_nuer1


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


       end subroutine get_f_nuer2


       end module trans_nuer2
