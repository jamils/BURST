
       module trans_nunubar


!------linkages.

!      used by main

!------remarks.

!      Module for neutrinos scattering on anti-neutrinos.


       use mainvar
       use transvar
       use interp


       implicit none


       contains


!-----------------------------------------------------------

       !subroutine scattering_nunubar(nutrans,eps_bins,p3_bins,nets,frss,neta,frsa)
       subroutine scattering_nunubar(nutrans,p3_bins,nets,frss,neta,frsa)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] none
!                  [function] none


!-----remarks.

!     Calculates scattering rates for \nu scattering on \bar{\nu} and
!     annihilation rates for \nu_i + \bar{\nu}_i -> \nu_j + \bar{\nu}_j.


!------modules.

       use gauss


       implicit none


       save


!------throughput variables.

       type (nuvar), intent(in), dimension(:,:) :: nutrans !nu occupation fractions and info.
       !type (bin_scheme), intent(in) :: eps_bins !epsilon abscissas and weights
       type (bin_scheme), intent(in) :: p3_bins !abscissas and weights for p3 integral
       real(dl), intent(out), dimension(:,:,:,:,:) :: nets !net scattering rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:,:) :: frss !forward-reverse summed scattering rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:) :: neta !net annihilation rates/tnmev**5
       real(dl), intent(out), dimension(:,:,:,:) :: frsa !forward-reverse summed annihilation rates/tnmev**5


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
       real(dl) weight3

       real(dl) forfact, revfact

       real(dl), dimension(3,2,3,2) :: nets1
       real(dl), dimension(3,2,3,2) :: frss1
       real(dl), dimension(3,2,3,2) :: nets2
       real(dl), dimension(3,2,3,2) :: frss2

       real(dl), dimension(3,3,2) :: neta1
       real(dl), dimension(3,3,2) :: frsa1
       real(dl), dimension(3,3,2) :: neta2
       real(dl), dimension(3,3,2) :: frsa2

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

       do n=1,2
         do m=1,3 !don't use sterile occupation fractions
           woccfrac(m,n,:) = nutrans(m,n)%occfrac(:)
         end do
       end do
       woccarray = 0._dl


       prefac = 128._dl*gfhbar/16._dl/(2._dl*pi)**3
       !testsum = 0._dl
       !testsum2 = 0._dl
       !lowsum = 0._dl
       !highsum = 0._dl

       startind = 1
       if (eps_bins%abscissas(1).eq.0._dl) startind = 2

       nets = 0._dl
       frss = 0._dl
       neta = 0._dl
       frsa = 0._dl
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
           kinemfact1 = eps_bins%weights(j)*p1
           nets1 = 0._dl
           frss1 = 0._dl
           neta1 = 0._dl
           frsa1 = 0._dl
           !write (116,'(A,1pe12.5)') 'p2 = ', p2
           !if (.not.(((p1.eq.0._dl).and.(p2.eq.0._dl)))) then
           !tempint = 0._dl
           if ((p1.ne.0._dl).and.(p2.ne.0._dl)) then
             glegai = p3_bins%abscissas
             glegwi = p3_bins%weights
             call cscgqf(1,glegai,glegwi,0._dl,p2,successflag)
             do k=1,size(glegai)
               p3 = glegai(k)
               weight3 = calc_weight_nunubar_small(p1,p3)
               p4 = p1 + p2 - p3
               kinemfact2 = glegwi(k)*weight3*p3**3
               !call interp_occfrac(eps_bins%abscissas,dlog(woccfrac), &
               !                               woccfracddot,p3,logwoccfrac)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p3,logwoccfrac)
               woccarray(:,:,3) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,3) = 1._dl/(dexp(p3) + 1._dl)
               !call interp_occfrac(eps_bins%abscissas,dlog(woccfrac), &
               !                               woccfracddot,p4,logwoccfrac)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
               woccarray(:,:,4) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,4) = 1._dl/(dexp(p4) + 1._dl)
               !testsum(i,j) = testsum(i,j) + glegwi(k)*weight3*p3**3
               call get_f_nunubar(woccarray,nets2,frss2,neta2,frsa2)
               nets1 = nets1 + nets2*kinemfact2
               frss1 = frss1 + frss2*kinemfact2
               neta1 = neta1 + neta2*kinemfact2
               frsa1 = frsa1 + frsa2*kinemfact2
               !tempint = tempint + frss2(1,2,1,1)*kinemfact2
               !write (116,*) p3, weight3*p3**3
               !write (116,*) p3, weight3*p3**3*frss2(1,2,1,1)
             end do !k; 0<p3<p2 loop
             glegai = p3_bins%abscissas
             glegwi = p3_bins%weights
             call cscgqf(1,glegai,glegwi,p2,(p1 + p2),successflag)
             do k=1,size(glegai)
               p3 = glegai(k)
               weight3 = calc_weight_nunubar_large(p1,p2,p3)
               p4 = p1 + p2 - p3
               kinemfact2 = glegwi(k)*weight3*p3**3
               !call interp_occfrac(eps_bins%abscissas,dlog(woccfrac), &
               !                               woccfracddot,p3,logwoccfrac)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p3,logwoccfrac)
               woccarray(:,:,3) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,3) = 1._dl/(dexp(p3) + 1._dl)
               !call interp_occfrac(eps_bins%abscissas,dlog(woccfrac), &
               !                               woccfracddot,p4,logwoccfrac)
               call interp_occfrac(eps_bins%abscissas,dlog(woccfrac) &
                                              ,p4,logwoccfrac)
               woccarray(:,:,4) = dexp(logwoccfrac(:,:))
               !woccarray(:,:,4) = 1._dl/(dexp(p4) + 1._dl)
               !woccfrac_exact = 1._dl/(dexp(p3) + 1._dl)
               !write (112,*) p3, (-woccfrac_exact + woccfrac3(1,1))/woccfrac_exact!, logwoccfrac, woccfrac3(1,1)
               !testsum(i,j) = testsum(i,j) + glegwi(k)*weight3*p3**3
               !call get_f_nunubar(i,woccarray,net2,frs2)
               call get_f_nunubar(woccarray,nets2,frss2,neta2,frsa2)
               nets1 = nets1 + nets2*kinemfact2
               frss1 = frss1 + frss2*kinemfact2
               neta1 = neta1 + neta2*kinemfact2
               frsa1 = frsa1 + frsa2*kinemfact2
               !write (116,*) p3, weight3*p3**3
               !write (116,*) p3, weight3*p3**3*frss2(1,2,1,1)
               !tempint = tempint + frss2(1,2,1,1)*kinemfact2
             end do !k; p2<p3<(p1 + p2) loop
           end if !.not.(p1 = 0 .and. p2 = 0)
           nets(1:3,:,:,:,i) = nets(1:3,:,:,:,i) + nets1(:,:,:,:)*kinemfact1
           frss(1:3,:,:,:,i) = frss(1:3,:,:,:,i) + frss1(:,:,:,:)*kinemfact1
           neta(1:3,:,:,i) = neta(1:3,:,:,i) + neta1(:,:,:)*kinemfact1
           frsa(1:3,:,:,i) = frsa(1:3,:,:,i) + frsa1(:,:,:)*kinemfact1
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
       nets = prefac*nets
       frss = prefac*frss
       neta = prefac*neta
       frsa = prefac*frsa

       !write (*,*) lowsum, highsum


       !deallocate(glegai)
       !deallocate(glegwi)

       !do i=1,nbins
       !  write (114,*) eps_bins%abscissas(i),((testsum(i,j) - 8._dl/9._dl*eps_bins%abscissas(j)**3) &
       !                                      *9._dl/8._dl/eps_bins%abscissas(j)**3,j=1,nbins)
       !  !write (117,*) eps_bins%abscissas(i),eps_bins%weights(i),net(1,1,1,1,i),frs(1,1,1,1,i) &
       !  !              ,(net(1,1,1,1,i)/frs(1,1,1,1,i))
       !end do
       !write (117,*) (nets(1,2,1,1,i)/frss(1,2,1,1,i),i=1,nbins)


       !deallocate(testsum)


       return


       end subroutine scattering_nunubar


!-------------------------------------------------------------------     
       
       real(dl) function calc_weight_nunubar_small(p1,p3)


!------linkages.

!      called by - [program] main
!      calls     - none


!------remarks.

!      Calculates weight function for \nu \nu scattering.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: p1 !first epsilon quantity
       real(dl), intent(in) :: p3 !third epsilon quantity


!------local variables.

       real(dl) fval !function value
       real(dl) psmall, plarge !epsilon values


!------procedure.


       psmall = min(p1,p3)
       plarge = max(p1,p3)
       fval = 4._dl/15._dl/plarge**3*(psmall**2 - 5._dl*p1*p3 + 10._dl*plarge**2)         

       !fail-safe:
       if (fval.lt.0._dl) fval = 0._dl

       calc_weight_nunubar_small = fval


       return


       end function calc_weight_nunubar_small


!-------------------------------------------------------------------     
       
       real(dl) function calc_weight_nunubar_large(p1,p2,p3)


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
       real(dl) prefac, first_term, second_term


!------procedure.


       if (p3.ge.(p1 + p2)) then
         fval = 0._dl
       else if (p3.eq.0._dl) then
         fval = 0._dl
       else
         second_term = 10._dl*(p1 - p3)**2 + 15._dl*(p1 - p3)*p2 + 6._dl*p2**2
         if (p1.ge.p3) then
           prefac = 4._dl*p2**3/15._dl/p1**3/p3**3
           fval = prefac*second_term
         else
           prefac = 4._dl/15._dl/p1**3/p3**3
           first_term = (p1 - p3)**5
           second_term = p2**3*second_term
           fval = prefac*(first_term + second_term)
         end if
       end if

       !fail-safe:
       if (fval.lt.0._dl) fval = 0._dl

       calc_weight_nunubar_large = fval


       return


       end function calc_weight_nunubar_large


!-----------------------------------------------------------

       !subroutine get_f_nunubar(i,woccarray,net,frs)
       subroutine get_f_nunubar(woccarray,nets,frss,neta,frsa)


!------linkages.

!      called by - [subroutine] scattering_nunubar
!      calls     - none

!------remarks.

!      Subroutine to assign nutrans values.


       implicit none


!------throughput variables.

       !integer i !index for net and frs
       real(dl), dimension(:,:,:), intent(in) :: woccarray !Array of occupation probs.
       !real(dl), dimension(:,:,:,:,:), intent(out) :: net !net values
       !real(dl), dimension(:,:,:,:,:), intent(out) :: frs !frs values
       real(dl), dimension(:,:,:,:), intent(out) :: nets !net s values
       real(dl), dimension(:,:,:,:), intent(out) :: frss !frs s values
       real(dl), dimension(:,:,:), intent(out) :: neta !net a values
       real(dl), dimension(:,:,:), intent(out) :: frsa !frs a values


!------local variables.

       integer n,nn,m,mm !indicies
       real(dl) flafact, forfact, revfact


!------procedure.


!------calculated values of frs and net rates.

       do n=1,2
         nn = 2 - n/2
         do m=1,3
           do mm=1,3
             flafact = 1._dl
             if (m.ne.mm) flafact = 0.25_dl
             forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(mm,nn,2)) &
                       *woccarray(mm,nn,3)*woccarray(m,n,4)
             revfact = woccarray(m,n,1)*woccarray(mm,nn,2) &
                       *(1._dl-woccarray(mm,nn,3))*(1._dl-woccarray(m,n,4))
             frss(mm,nn,m,n) = flafact*(forfact + revfact)
             nets(mm,nn,m,n) = flafact*(forfact - revfact)
             if (m.ne.mm) then
               forfact = (1._dl-woccarray(m,n,1))*(1._dl-woccarray(m,nn,2)) &
                         *woccarray(mm,nn,3)*woccarray(mm,n,4)
               revfact = woccarray(m,n,1)*woccarray(m,nn,2) &
                         *(1._dl-woccarray(mm,nn,3))*(1._dl-woccarray(mm,n,4))
               frsa(mm,m,n) = 0.25_dl*(forfact + revfact)
               neta(mm,m,n) = 0.25_dl*(forfact - revfact)
             end if
           end do !mm
         end do !m
       end do !n


       return


       end subroutine get_f_nunubar


       end module trans_nunubar
