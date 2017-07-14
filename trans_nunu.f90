
       module trans_nunu


!------linkages.

!      used by trans_init and trans_evolve

!------remarks.

!      Module for neutrino scattering on neutrinos.


       use mainvar
       use transvar
       use interp


       implicit none


       contains


!-----------------------------------------------------------

       !subroutine scattering_nunu(nutrans,eps_bins,p3_bins,net,frs)
       subroutine scattering_nunu(nutrans,p3_bins,net,frs)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] none
!                  [function] none


!-----remarks.

!     Calculates scattering rates for \nu scattering on \nu and
!     \bar{nu} scattering on \bar{\nu}.


!------modules.

       use gauss


       implicit none


       save


!------throughput variables.

       type (nuvar), intent(in), dimension(:,:) :: nutrans !nu occupation fractions and info.
       !type (bin_scheme), intent(in) :: eps_bins !epsilon abscissas and weights
       type (bin_scheme), intent(in) :: p3_bins !abscissas and weights for p3 integral
       real(dl), intent(out), dimension(:,:,:,:,:) :: net !net scattering rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:,:) :: frs !forward-reverse summed scattering rates/tnmev**5


!------local variables.

       integer i,j,k,m,mm,n
       !integer nbins
       integer sumind

       real(dl) p1, p2, p3, p4
       real(dl) prefac
       real(dl) kinemfact1, kinemfact2

       !real(dl) dlogwoccfrac1
       !real(dl) dlogwoccfrac2
       real(dl), dimension(3,2) :: logwoccfrac

       real(dl), dimension(3,2,4) :: woccarray !array for calculating F

       logical successflag

       !real(dl), dimension(:), allocatable :: glegai
       !real(dl), dimension(:), allocatable :: glegwi

       real(dl) weightinterp
       real(dl) weight3, psmall, plarge

       real(dl) forfact, revfact

       real(dl), dimension(3,2,3,2) :: net1
       real(dl), dimension(3,2,3,2) :: frs1
       real(dl), dimension(3,2,3,2) :: net2
       real(dl), dimension(3,2,3,2) :: frs2

       !real(dl) tempint

       integer startind !starting index


!------Testing quantities.

       !real(dl) lowsum, highsum
       !real(dl), dimension(:,:), allocatable :: testsum
       !real(dl) testsum2, woccfrac_exact


!------procedure.


       !allocate(glegai(size(p3_bins%abscissas))) !assign dimension of glegai
       !allocate(glegwi(size(p3_bins%abscissas))) !assign dimension of glegwi
       

       !nbins = size(nutrans(1,1)%occfrac)

       !allocate(testsum(nbins,nbins)) !assign dimension of testsum

       !do i=1,nbins
       !  do n=1,2
       !    do m=1,3 !don't use sterile occupation fractions
       !      woccfrac(m,n,i) = nutrans(m,n)%occfrac(i)
       !    end do
       !  end do
       !end do
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

       prefac = 128._dl*gfhbar/32._dl/(2._dl*pi)**3
       !testsum = 0._dl
       !testsum2 = 0._dl
       !lowsum = 0._dl
       !highsum = 0._dl

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2


       net = 0._dl
       frs = 0._dl
       sumind = 1
       do i=startind,nbins
       !do i=2,nbins,10
       !do i=2,2
         !p1 = eps_bins%abscissas(i)
         !write (118,'(A,1pe12.5)') 'p1 = ', p1
         do j=startind,nbins
         !do j=2,nbins,10
         !do j=nbins,nbins
#ifdef prllel
         if (rank.eq.mod(sumind,nprocs)) then
#endif
           p1 = eps_bins%abscissas(i)
           p2 = eps_bins%abscissas(j)
           woccarray(:,:,1) = woccfrac(:,:,i)
           woccarray(:,:,2) = woccfrac(:,:,j)
           kinemfact1 = eps_bins%weights(j)*p1*p2**3
           net1 = 0._dl
           frs1 = 0._dl
           if (i.ge.j) then
             psmall = p2
             plarge = p1
           else if (j.gt.i) then
             psmall = p1
             plarge = p2
           end if
           !write (116,'(A,1pe12.5)') 'p2 = ', p2
           !if (.not.(((p1.eq.0._dl).and.(p2.eq.0._dl)))) then
           !tempint = 0._dl
           if ((p1.ne.0._dl).and.(p2.ne.0._dl)) then
             glegai = p3_bins%abscissas
             glegwi = p3_bins%weights
             call cscgqf(1,glegai,glegwi,0._dl,psmall,successflag)
             do k=1,size(glegai)
               p3 = glegai(k)
               weight3 = calc_weight_nunu_ends(p1,p2,p3)
               p4 = p1 + p2 - p3
               kinemfact2 = glegwi(k)*weight3
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p3,logwoccfrac)
               woccarray(:,:,3) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,3) = 1._dl/(dexp(p3) + 1._dl)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
               woccarray(:,:,4) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,4) = 1._dl/(dexp(p4) + 1._dl)
               !testsum(i,j) = testsum(i,j) + glegwi(k)*weight3
               call get_f_nunu(woccarray,net2,frs2)
               net1 = net1 + net2*kinemfact2
               frs1 = frs1 + frs2*kinemfact2
               !tempint = tempint + frs2(1,1,1,1)*kinemfact2
               !write (116,*) p3, weight3*frs2(1,1,1,1)
             end do !k; 0<p3<psmall loop
           end if !.not.(p1 = 0 .and. p2 = 0)
           if (i.ne.j) then
             glegai = p3_bins%abscissas
             glegwi = p3_bins%weights
             weight3 = calc_weight_nunu_mid(p1,p2)
             call cscgqf(1,glegai,glegwi,psmall,plarge,successflag)
             do k=1,size(glegai)
               p3 = glegai(k)
               p4 = p1 + p2 - p3
               kinemfact2 = glegwi(k)*weight3
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p3,logwoccfrac)
               woccarray(:,:,3) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,3) = 1._dl/(dexp(p3) + 1._dl)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
               woccarray(:,:,4) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,4) = 1._dl/(dexp(p4) + 1._dl)
               !woccfrac_exact = 1._dl/(dexp(p3) + 1._dl)
               !testsum(i,j) = testsum(i,j) + glegwi(k)*weight3
               !call get_f_nunu(i,woccarray,net2,frs2)
               call get_f_nunu(woccarray,net2,frs2)
               net1 = net1 + net2*kinemfact2
               frs1 = frs1 + frs2*kinemfact2
               !write (116,*) p3, weight3*frs2(1,1,1,1)
               !tempint = tempint + frs2(1,1,1,1)*kinemfact2
             end do !k; pmall<p3<plarge loop
           end if !(i.ne.j)
           !if (.not.(((p1.eq.0._dl).and.(p2.eq.0._dl)))) then
           if ((p1.ne.0._dl).and.(p2.ne.0._dl)) then
             glegai = p3_bins%abscissas
             glegwi = p3_bins%weights
             call cscgqf(1,glegai,glegwi,plarge,(p1 + p2),successflag)
             do k=1,size(glegai)
               p3 = glegai(k)
               weight3 = calc_weight_nunu_ends(p1,p2,p3)
               p4 = p1 + p2 - p3
               kinemfact2 = glegwi(k)*weight3
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p3,logwoccfrac)
               woccarray(:,:,3) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,3) = 1._dl/(dexp(p3) + 1._dl)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
               woccarray(:,:,4) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,4) = 1._dl/(dexp(p4) + 1._dl)
               !testsum(i,j) = testsum(i,j) + glegwi(k)*weight3
               call get_f_nunu(woccarray,net2,frs2)
               net1 = net1 + net2*kinemfact2
               frs1 = frs1 + frs2*kinemfact2
               !write (116,*) p3, weight3*frs2(1,1,1,1)
               !tempint = tempint + frs2(1,1,1,1)*kinemfact2
             end do !k; plarge<p3<p1+p2 loop
           end if
           net(1:3,:,:,:,i) = net(1:3,:,:,:,i) + net1(:,:,:,:)*kinemfact1
           frs(1:3,:,:,:,i) = frs(1:3,:,:,:,i) + frs1(:,:,:,:)*kinemfact1
           !if ((i.eq.100).and.(j.eq.99)) write (117,*) net(1,1,1,1,i), frs(1,1,1,1,i)
           !write (116,*)
#ifdef prllel
         end if
#endif
           sumind = sumind + 1
           !tempint = tempint*p2**3
           !write (118,*) p2, tempint
         end do !j; p_2 loop
         !write (118,*)
       end do !i; p_1 loop
       net = prefac*net
       frs = prefac*frs

       !write (*,*) lowsum, highsum


       !deallocate(glegai)
       !deallocate(glegwi)

       !do i=1,nbins
       !  write (114,*) eps_bins%abscissas(i),((testsum(i,j) - 8._dl/3._dl)*3._dl/8._dl,j=1,nbins)
       !  !write (117,*) eps_bins%abscissas(i),eps_bins%weights(i),net(1,1,1,1,i),frs(1,1,1,1,i) &
       !  !              ,(net(1,1,1,1,i)/frs(1,1,1,1,i))
       !end do
       !write (117,*) (net(1,1,1,1,i)/frs(1,1,1,1,i),i=1,nbins)


       !deallocate(testsum)


       return


       end subroutine scattering_nunu


!-------------------------------------------------------------------     
       
       real(dl) function calc_weight_nunu_ends(p1,p2,p3)


!------linkages.

!      called by - [program] main
!      calls     - none


!------remarks.

!      Calculates weight function for \nu \nu scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: p2 !second epsilon quantity
       real(dl), intent(in) :: p3 !third epsilon quantity


!------local variables.

       real(dl) fval !function value
       real(dl) prefac, first_term, second_term !terms


!------procedure.


       if (p3.ge.(p1 + p2)) then
         fval = 0._dl
       else if (p3.eq.0._dl) then
         fval = 0._dl
       else
         first_term = (p1 + p2)**4
         if ((p1+p2).ge.(2._dl*p3)) then
           prefac = 4._dl*p3/15._dl/p1**3/p2**3
           second_term = -(p1 + p2 - p3)*(p1 + p2 - 2._dl*p3)*((p1 + p2)**2 + 3._dl*p3*(p1 + p2 - p3))
         else
           prefac = 4._dl*(p1 + p2 - p3)/15._dl/p1**3/p2**3
           second_term = (p1 + p2 - 2._dl*p3)*p3*((p1 + p2)**2 + 3._dl*p3*(p1 + p2 - p3))
         end if
         fval = prefac*(first_term + second_term)
       end if

       !fail-safe:
       if (fval.lt.0._dl) fval = 0._dl

       calc_weight_nunu_ends = fval


       return


       end function calc_weight_nunu_ends


!-------------------------------------------------------------------     
       
       real(dl) function calc_weight_nunu_mid(p1,p2)


!------linkages.

!      called by - [program] main
!      calls     - none


!------remarks.

!      Calculates weight function for \nu \nu scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: p2 !second epsilon quantity


!------local variables.

       real(dl) fval !function value
       real(dl) psmall, plarge !small and large values of (p1,p2)


!------procedure.


       psmall = min(p1,p2)
       plarge = max(p1,p2)
       fval = 4._dl/15._dl/plarge**3*(psmall**2 + 5._dl*p1*p2 + 10._dl*plarge**2)         

       !fail-safe:
       if (fval.lt.0._dl) fval = 0._dl

       calc_weight_nunu_mid = fval


       return


       end function calc_weight_nunu_mid


!-----------------------------------------------------------

       subroutine get_f_nunu(woccarray,net,frs)


!------linkages.

!      called by - [program] scattering_nunu
!      calls     - none

!------remarks.

!      Subroutine to assign nutrans values.


       implicit none


!------throughput variables.

       real(dl), dimension(:,:,:), intent(in) :: woccarray !Array of occupation probs.
       real(dl), dimension(:,:,:,:), intent(out) :: net !net values
       real(dl), dimension(:,:,:,:), intent(out) :: frs !frs values


!------local variables.

       integer n,m,mm !indicies
       real(dl) flafact, forfact, revfact


!------procedure.


!------calculated values of frs and net rates.

       do n=1,2
         do m=1,3
           do mm=1,3
             flafact = 1._dl
             if (m.ne.mm) flafact = 0.5_dl
             forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(mm,n,2)) &
                       *woccarray(m,n,3)*woccarray(mm,n,4)
             revfact = woccarray(m,n,1)*woccarray(mm,n,2) &
                       *(1._dl-woccarray(m,n,3))*(1._dl-woccarray(mm,n,4))
             frs(mm,n,m,n) = flafact*(forfact + revfact)
             net(mm,n,m,n) = flafact*(forfact - revfact)
           end do !mm
         end do !m
       end do !n


       return


       end subroutine get_f_nunu


       end module trans_nunu
