
       module sdratio_length


!------linkages.

!      used by - [subroutine] get_neff


!------remarks.

!      Contains code to do integrals for r_s and r_d


!------modules.


       use sdratiovar


       implicit none


!------Boltzmann eqn. variables.

       real(dl), dimension(:), allocatable :: xec !total free-electron fraction
       !real(dl), dimension(:,:), allocatable, target :: xei !array of free-electron fractions
       real(dl), dimension(:), allocatable :: ascarray !array of scale factors for xec
       !real(dl), dimension(:), allocatable, target :: tbar !baryon temperature

       real(dl), dimension(:,:), allocatable, target :: ratioy !pointer to evolve xei and tbar
       real(dl), dimension(:,:), pointer :: xei !array of free-electron fractions
       real(dl), dimension(:), pointer :: tbar !baryon temperature


       real(dl) adcod
       !integer :: counter = 1


       contains


!--------------------------------------------------------------------

       subroutine sdratio_length_calc(thryflag,difflenc,ratioarrayout,soundhor)


!------linkages.

!      called by - [subroutine] sdratio
!      calls     - [subroutine] get_sel_recom_history, find_adc
!                  [function] occfraclagint


!------remarks.

!      Performs integration over a to retrieve sound horizon and diffusion length.
!      If thryflag = .true. then soundhor is not needed.


!------modules.


       use xe_history
       use opt_depth_mod
       use ang_diam_mod
       use interp


       implicit none


       save


!------throughput variables.

       !real(dl), intent(in) :: aend !upper limit of integral
       logical, intent(in) :: thryflag
       real(dl), intent(out) :: difflenc
       real(dl), intent(out), dimension(:) :: ratioarrayout
       real(dl), optional, intent(out) :: soundhor


!------local variables.

       real(dl) asc !integration variable, scale factor
       real(dl) dasc !integration variable width
       real(dl) apri !pervious value of scale factor
       real(dl) astep !step size for RK4 integration

       real(dl) rhobint !baryon density at asc
       real(dl) rhogint !photon density at asc
       real(dl) rint ! = 3*rhobint/4/4rhogint

       real(dl) ggs !energy statistical weight at a = 0
       real(dl) tdil !Dilution factors for neutrino temperatures
       real(dl) hubrate !Hubble rate

       real(dl) intfacts !integrand for sound horizon
       real(dl) intfactdc !integrand for square of the diffusion length
       real(dl) difflenc2 !square of the diffusion length

       integer m,n,i,j,k !indicies
       integer brnum !Used for brcoeffs in integration
       integer intnum2 !number of steps in integral

       real(dl) tempmev !Photon temperature in MeV
       real(dl) nh !number density of H ions/atoms
       real(dl) nhe !number density of He ions/atoms
       real(dl) rhs !right hand side of Saha eqn.

       real(dl) xetot !Total ionization fraction for H and He (single and double)
       real(dl) difflenctemp, soundhortemp, ratiosdtemp
       real(dl) lagint

       real(dl) ain, tau_od
       real(dl) alow, ahigh
       integer counter
       real(dl) d_a

       real(dl), dimension(:), allocatable :: xecddot !Second derivatives of xec
       real(dl), dimension(4) :: xeistart !starting array of values for xei
       real(dl), dimension(3) :: xeitemp !array of values for xei at single value of a
       real(dl), dimension(1001) :: xectemp
       real(dl) aguess


!------procedure.

       !if (.not.thryflag) adc = 9.1923934718086964e-04_dl !for 0.23 eV mass
       !if (.not.thryflag) adc = 9.1898085634656138e-04_dl !for 1.0 eV mass

       xec = 0._dl !total electron fraction for boltzmann method w/ corr.
       xei = 0._dl !array of electron fractions for boltzmann method w/ corr.
       tbar = 0._dl !baryon temperature

       !adcod = 9.1920328474044796E-004_dl !only used for testing purposes.

       if (.not.thryflag) then
         xeistart(1) = xeimax(1)
         xeistart(2) = 0._dl
         xeistart(3) = xeimax(3)
         xeistart(4) = 0._dl
         call get_sel_recom_history(ratioy,xeistart,ascarray,0._dl,1._dl,.false.)
         !write (106,*) 'Y_P = ', ypxe, ' \omega_b = ', obh2
         !write (106,*) lastind
         !do i=1,lastind
         !  write (105,*) ascarray(i), t0mev/ascarray(i), ratioy(i,4)
         !  write (106,*) ascarray(i), sum(ratioy(i,1:3)*xeiz(1:3))
         !  write (106,*) ascarray(i), (ratioy(i,j),j=1,3)
         !end do
         !write (106,'(/)')
         !xeistart(1) = xeimax(1)
         !xeistart(2) = 0.99_dl*xeimax(2)
         !xeistart(2) = 0._dl
         !xeistart(3) = 0._dl
         !xeistart(3) = xeimax(3)
         !xeistart(4) = t0mev/atranshe21
         !xeistart(4) = 0._dl
         !call get_sel_recom_history(xei,xeistart,ascarray,atranshe21,1._dl,.false.)
         !call get_sel_recom_history(ratioy,xeistart,ascarray,atranshe21,1._dl,.false.)
         !call get_sel_recom_history(ratioy,xeistart,ascarray,0._dl,adcod,.true.)
       else
         xeistart(1) = xeimax(1)
         xeistart(2) = 0._dl
         xeistart(3) = xeimax(3)
         xeistart(4) = 0._dl
         !call get_sel_recom_history(xei,xeistart,ascarray,0._dl,adcod,.true.)
         call get_sel_recom_history(ratioy,xeistart,ascarray,0._dl,adcod,.true.)
       end if

       do m=1,3
         xec(:) = xec(:) + xei(:,m)*xeiz(m)
       end do
       !if (neffval.eq.3._dl) xectemp = xec(1:1001)

       if (.not.thryflag) then
         !aguess = aend
         call find_adc(ascarray,xec,hubparm,adcod)
         !adcod = 9.1920328474044796E-004_dl !only used for testing purposes.
         !adcod = adcpl
         xeistart(1) = xeimax(1)
         xeistart(2) = 0._dl
         xeistart(3) = xeimax(3)
         xeistart(4) = 0._dl
         xei = 0._dl
         tbar = 0._dl !baryon temperature
         xec = 0._dl
         call get_sel_recom_history(ratioy,xeistart,ascarray,0._dl,adcod,.true.)
         do m=1,3
           xec(:) = xec(:) + xei(:,m)*xeiz(m)
         end do
         !do i=1,intnum
         !  write (101,*) ascarray(i), xec(i)!, t0mev/ascarray(i), tbar(i)
         !end do
       end if
       !d_a = ang_diam(adcod,hubparm)
       !write (*,*) 'Angular Diameter: ', d_a

       !if ((neffval.eq.3._dl).or.(.not.thryflag)) then
       !if (.not.thryflag) then
       !  !write (*,*) 'Planck'
       !  !write (*,*) 'a_{dc} = ', aend
       !  !tau_od = opt_depth(aend,ascarray,xec,hubparm)
       !  !write (*,*) 'Optical Depth: ', tau_od
       !  !d_a = ang_diam(aend,hubparm)
       !  !write (*,*) 'Angular Diameter: ', d_a

       !  !do i=1,lastind
       !  !  write (101,*) ascarray(i), xec(i)
       !  !end do
       !  !call flush(101)
       !  alow = 0.9_dl*aend
       !  do
       !    tau_od = opt_depth(alow,ascarray,xec,hubparm)
       !    if (tau_od.gt.1._dl) then
       !      exit
       !    else
       !      alow = 0.9_dl*alow
       !    end if
       !  end do
       !  ahigh = 1.1_dl*aend
       !  do
       !    tau_od = opt_depth(ahigh,ascarray,xec,hubparm)
       !    if (tau_od.lt.1._dl) then
       !      exit
       !    else
       !      ahigh = 1.1_dl*ahigh
       !    end if
       !  end do
       !  ain = aend
       !  counter = 1
       !  do
       !    tau_od = opt_depth(ain,ascarray,xec,hubparm)
       !    if ((abs(tau_od-1._dl).lt.1.e-06_dl).or.(counter.eq.100)) then
       !      exit
       !    else
       !      if (tau_od.gt.1._dl) then
       !        alow = ain
       !        ain = (ain + ahigh)/2._dl
       !      else
       !        ahigh = ain
       !        ain = (ain + alow)/2._dl
       !      end if
       !    end if
       !    counter = counter + 1
       !  end do
       !  write (*,*) 'a_{dc} = ', ain
       !  write (*,*) 'Optical Depth:', tau_od
       !  write (*,*) 'counter:', counter
       !  d_a = ang_diam(ain,hubparm)
       !  write (*,*) 'Angular Diameter:', d_a
       !end if

       !no contribution from first element of integral:
       if (present(soundhor)) then
         soundhor = 0._dl !numerically integrated sound horizon
         ggs = 1._dl
         do n=1,2
           do m=1,size(hubparm%nuparm, 1)
             !19Apr2015 EG: dilfact not needed as it's built into nuvar%occfrac:
             !tdil = hubparm%nuparm(m,n)%dilfact
             !ggs = ggs + tdil**4*(4._dl/11._dl)**(4._dl/3._dl)*7._dl/8._dl/2._dl
             !15Apr2015 EG: Change occfraclagint to different function.
             !lagint = occfraclagint(hubparm%nuparm(m,n)%occfrac,0._dl)
             !ggs = ggs + tdil**4*(4._dl/11._dl)**(4._dl/3._dl)/2._dl*lagint*15._dl/pi**4
             lagint = occfracint(hubparm%nuparm(m,n)%occfrac,0._dl)
             !ggs = ggs + (4._dl/11._dl)**(4._dl/3._dl)/2._dl*lagint*15._dl/pi**4
             ggs = ggs + hubparm%nuparm(m,n)%dilfact**4/2._dl*lagint*15._dl/pi**4
           end do
         end do
         ggs = ggs + hubparm%orad*hubparm%rhoc0/t0mev**4*15._dl/pi**2
       end if

       !no contribution from first element of integral:
       !numerically integrated square of diffusion length...
       !(boltzmann-solver with correction method):
       difflenc2 = 0._dl

       !width of integration variable:
       dasc = adcod/real(intnum-1, kind=dl)

       !if (.not.thryflag) then
       !  allocate(xecddot(lastind)) !Declare size of xecddot
       !  call cubic_spline(ascarray(1:lastind),xec(1:lastind),xecddot,0._dl,0._dl)
       !end if


       do i=1,intnum

         brnum = mod(i-2,4) + 1
         !indexed scale factor a:
         !asc = real(i-1, kind=dl)/real(intnum-1, kind=dl)*aend
         asc = real(i-1, kind=dl)/real(intnum-1, kind=dl)*adcod
         !baryon energy density at a = asc:
         !rhobint = rhob0*(a0/asc)**3
         !photon energy density at a = asc:
         !rhogint = rhog0*(a0/asc)**4
         !calculate hubrate at a = asc:
         if (asc.ne.0._dl) then
           hubrate = hubcalc(asc,hubparm)
         else
           hubrate = 0._dl
         end if
         !r at a = asc:
         !rint = 0.75_dl*rhobint/rhogint
         rint = 0.75_dl*rhob0/rhog0*asc/a0

         !integrate the sound horizon:
         if (present(soundhor)) then
           if (i.eq.1) then
             intfacts = brends(1)*dasc/sqrt(8._dl*pi)*mpl/sqrt(rhog0*ggs)
           else if (i.eq.intnum) then
             intfacts = dasc/asc**2/hubrate/sqrt(3._dl*(1._dl + rint))
             intfacts = brends(2)*intfacts
           else if (i.gt.intnum) then
             intfacts = 0._dl
           else
             intfacts = dasc/asc**2/hubrate/sqrt(3._dl*(1._dl + rint))
             intfacts = brcoeffs(brnum)*intfacts
           end if
           soundhor = soundhor + intfacts !in (Mpc*MeV)^-1 (comoving)
         end if

         !calcualate total free-electron fraction:
         !if (((neffval.eq.3._dl).or.(.not.thryflag)).and.(xeflag)) then
         !if ((.not.thryflag).and.(xeflag)) then
         !  !Going to need to grab identical recom history until ~10^-3
         !  if (i.eq.1) then
         !    xetot = 1._dl
         !  else if (asc.lt.ascarray(1)) then
         !    call saha_net(asc,xeitemp)
         !    xetot = sum(xeitemp*xeiz)
         !  else if (asc.gt.ascarray(lastind)) then !should never go into this else statement
         !    !xetot = xec(size(xec))
         !    xetot = xec(lastind)
         !  else
         !    !call cubic_splint(ascarray,xec,xecddot,asc,xetot)
         !    call cubic_splint(ascarray(1:lastind),xec(1:lastind),xecddot,asc,xetot)
         !  end if
         !else
         !  xetot = xec(i)
         !end if
         xetot = xec(i)
         !integrate the diffusion length squared (boltzmann-corr.):
         if (hubrate.ne.0._dl) then
           intfactdc = pi**2*dasc/6._dl/ne0/xetot/sigmat/hubrate/ &
                       (1._dl + rint)*(rint**2/(1._dl + rint) + 16._dl/15._dl)
                       !in (Mpc*MeV)^-2 (comoving)
         else
           intfactdc = 0._dl
         end if
         !if ((neffval.eq.3._dl).and.(xeflag)) then
         !  write (97,*) asc, xetot, xei(i,3)
         !end if
         !if ((.not.thryflag).and.(xeflag)) then
         !  write (97,*) asc, xetot, xectemp(i), (xetot - xectemp(i))/xectemp(i)
         !end if

         if (i.eq.1) then
           intfactdc = 0._dl
         else if (i.eq.intnum) then
           intfactdc = brends(2)*intfactdc
         else if (i.gt.intnum) then
           intfactdc = 0._dl
         else
           intfactdc = brcoeffs(brnum)*intfactdc
         end if
         difflenc2 = difflenc2 + intfactdc

         if ((i.ne.1).and.(brnum.eq.4).and.(neffevflag)) then
           if (i.ne.intnum) then
             difflenctemp = sqrt(difflenc2 + intfactdc*(brends(2)/brcoeffs(brnum) - 1._dl))*convfact1
             if ((thryflag)) then
               soundhortemp = shcalc(asc,neffval)
             else
               soundhortemp = (soundhor + intfacts*(brends(2)/brcoeffs(brnum) - 1._dl))*convfact1
             end if
             ratioarrayout((i-1)/4) = soundhortemp/difflenctemp
             if ((neffval.eq.3._dl).or.(.not.thryflag)) then
               write (100,*) asc, soundhortemp, difflenctemp
             end if
           end if
         end if

       end do

       !convert to comoving dimensionless units:
       if (present(soundhor)) then
         soundhor = soundhor*convfact1
       end if
       difflenc = sqrt(difflenc2)*convfact1


       !if ((neffval.eq.3._dl).and.(xeflag)) then
       !  write (97,*) neffval, thryflag
       !  write (97,'(/)') 
       !  write (98,'(/)')
       !end if
       !if ((.not.thryflag).and.(xeflag)) then
       !  write (97,*) thryflag
       !  write (97,'(/)') 
       !  write (98,'(/)')
       !end if


       !if (.not.thryflag) deallocate(xecddot)


       end subroutine sdratio_length_calc


!--------------------------------------------------------------------

       real(dl) function nufdm(x,ms,deg)
       

!------linkages.

!      called by - [function] rhonulagint
!      calls     - [function] none


!------remarks.

!      Analytic form of expression for neutrino-energy distribution.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !integration variable
       real(dl), intent(in) :: ms !nu mass/nu temp.
       real(dl), intent(in) :: deg !nu degneracy parameter


!------procedure.

       !nufdm = sqrt(x**2 + ms**2)*x**2/(exp(x - deg) + 1._dl) !For 15g rule
       !nufdm = x**2*sqrt(x**2 + ms**2)*exp(x)/(exp(x - deg) + 1._dl) !For 64e rule
       nufdm = sqrt(x**2 + ms**2)/x*dexp(deg)/(1._dl + dexp(-(x - deg))) !For 64 rule
       !nufdm = x**2*sqrt(x**2 + ms**2)/(exp(x - deg) + 1._dl) !For rectangle rule
       
       return


       end function nufdm


       end module sdratio_length
