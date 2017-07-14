
       subroutine sdratio(obh2in,neffparm,ypri,tneff,nuvarin)
       
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] sdratio_launch, sdratio_length_calc
!                  [function] rhonulagint


!-----remarks.

!     Computes the sound horizon and diffusion length at surface of last scattering,
!     relating the two quantities to Neff


!------modules.

       use sdratiovar
       use sdratio_length
       use interp
       use massvar
       use xe_history
       !use ang_diam_mod


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: obh2in !input baryon number
       real(dl), intent(in) :: neffparm !parameter for Neff
       real(dl), intent(in) :: ypri !primordial 4He mass fraction
       real(dl), intent(inout) :: tneff !Neff target in; Neff calculated out
       type (nuvar), intent(inout), dimension(:,:), target :: nuvarin !nu occupation fractions and info.


!------local variables.

       real(dl) aeq !epoch of matter-radn. equality

       !real(dl) omegag !contribution from photons
       real(dl) omeganu !contribution from nu
       !real(dl) omegarad !contribution from extra relativistic energy
       real(dl) omegal !contribution from dark energy

       !real(dl), dimension(:), allocatable :: equilfrac !equilibrium occupation fractions.

       !integer neffnum !array size
       !real(dl) neffwidth !width of Neff range

       real(dl) soundhor
       real(dl) difflenc
       real(dl) ratiosd
       !real(dl), dimension(:,:), allocatable :: ratioarray
       !real(dl), dimension(:), allocatable :: neffarray
       !!!real(dl), dimension(:,:), allocatable :: ratioddot
       real(dl) tneffpred

       real(dl) kdc

       real(dl) mirrorwidth !width in Neff range around endpoints
       real(dl), dimension(4) :: neffend !end point Neff values
       !real(dl), dimension(:,:), allocatable :: ratioend !end point ratio values
       !real(dl), dimension(:), allocatable :: ratiocalc !calculated ratio values for input cosmology
       !real(dl), dimension(:), allocatable :: dratio !difference in ratio values
       !real(dl), dimension(:), allocatable :: dneff1 !Derivative of neff wrt ratio at left endpoint
       !real(dl), dimension(:), allocatable :: dneff2 !Deritative of neff wrt ratio at right endpoint

       real(dl) tempnuev
       real(dl) ms !nu mass/temp.
       real(dl) rhonui0

       real(dl) ratiotest
       real(dl) soundhortest
       real(dl) difflenctest

       integer i,j,nind,m,n
       logical thryflag
       !integer neffevnum !dimension for number of points to store when neffevflag = .true.

       !real(dl) summnulow !lowest value of \Sigma m_\nu
       !real(dl) summnuhigh !highest value of \Sigma m_\nu

       integer summnuind
       real(dl) summnu !\Sigma m_\nu
       !real(dl) masstol !tolerance to accept massl
       !logical hierflag !Flag for neutrino mass hierarchy

       real(dl), dimension(3) :: acmass !active nu masses
       real(dl) deltaneff

       real(dl) d_a !Angular Diameter distance to last scattering
       !real(dl) adcpl
       real(dl) ratiosdcalc


!------data file.


10     format (6x,'neff',8x,'kdc',9x,'rsi',9x,'rsx',8x,'ratio',/,60('-'))

25     format (1p10e12.3)



!------procedure.

       !write (*,*) nuvarin(1,1)%occfrac(1:3)

       call sdratio_launch(obh2in,neffparm,ypri,tneff)



       do summnuind=1,summnuprec

!------calculate ratio based on input cosmology.

       thryflag = .false.
       neffvalxe = 5._dl

       !write (*,*) neffparm, drhorad0, tdil, orad

       hubparm%onflag = .false.
       allocate(hubparm%nuparm(size(nuvarin,1),2)) !assign dimension of huparm%nuparm with non-standard cosmology


       !Determine \Sigma m_\nu:
       if (summnuprec.eq.1) then
         summnu = summnulow
       else
         summnu = real(summnuind-1, kind=dl)*(summnuhigh - summnulow)/real(summnuprec - 1, kind=dl) &
                  + summnulow
       end if

#ifdef prllel
       if (rank.eq.0) then
#endif

       !write (92,*) summnu

#ifdef prllel
       end if
#endif

       !active neutrino masses in eV:
       call get_numass(summnu,hierflag,acmass,masstol)
       !deltaneff = 0._dl
       !do m=1,3
       !  deltaneff = deltaneff + 5._dl/7._dl/pi**2*(11._dl/4._dl)**(2._dl/3._dl) &
       !                          *(acmass(m)/temp0ev)**2*(adcpl/a0)**2
       !end do
       !write (*,*) summnu, deltaneff
       do m=1,3
         nuvarin(m,:)%mass = acmass(m)
       end do

       !nuvarin(4,:)%mass = stmass
       hubparm%nuparm => nuvarin


       !contribution from neutrinos at current epoch:
       omeganu = 0._dl
       do n=1,2
         do m=1,size(hubparm%nuparm, 1)
           !tempnuev = temp0ev*(4._dl/11._dl)**(1._dl/3._dl)
           tempnuev = temp0ev*hubparm%nuparm(m,n)%dilfact
           ms = hubparm%nuparm(m,n)%mass/tempnuev
           rhonui0 = 1._dl/2._dl/pi**2*tempnuev**4 &
                     !*1.e-24_dl*occfraclagint(hubparm%nuparm(m,n)%occfrac, ms)
                     *1.e-24_dl*occfracint(hubparm%nuparm(m,n)%occfrac, ms)
           omeganu = omeganu + rhonui0/hubparm%rhoc0
         end do
       end do

       omegal = 1._dl - hubparm%om - hubparm%og - omeganu - hubparm%orad

       !energy density of dark energy at current epoch in units of MeV**4:
       rhode = omegal*hubparm%rhoc0
       hubparm%ol = omegal

       call sdratio_length_calc(thryflag,difflenc,ratiocalc,soundhor)       

       !d_a = ang_diam(adcpl,hubparm)
       !write (*,*) 'Angular Diameter:', d_a

       ratiosd = soundhor/difflenc
       ratiocalc(size(ratiocalc)) = ratiosd
       ratiosdcalc = ratiosd

#ifdef prllel
       if (rank.eq.0) then
#endif
       !write (*,*) soundhortest, difflenctest, ratiotest
       !write (*,*) soundhor, difflenc, ratiosd
       !write (*,*) soundhortest,soundhor
       !write (*,*) difflenctest,difflenc,pi/difflenc
#ifdef prllel
       end if !(rank.eq.0)
#endif


!------calculate theory points.

       ratioarray = 0._dl
       neffarray = 0._dl
       hubparm%orad = 0._dl
       thryflag = .true.
       hubparm%onflag = .true.

       do nind=1,neffnum !end do at 100

         neffval = real(nind-1, kind=dl)*(neffrangeh - neffrangel)/real(neffnum - 1, kind=dl) &
                   + neffrangel

         neffvalxe = neffval

         !contribution from purely relativistic neutrinos:
         !omeganu = omegag*7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval
         omeganu = hubparm%og*7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval
         hubparm%on = omeganu
         !omegal = 1._dl - hubparm%om - omegag - omeganu !contribution from dark energy
         omegal = 1._dl - hubparm%om - hubparm%og - omeganu !contribution from dark energy
         !energy density of dark energy at current epoch in units of MeV**4:
         rhode = omegal*hubparm%rhoc0
         hubparm%ol = omegal
         !if (nind.eq.6) then
         !  write (*,*) neffval
         !  write (*,*) obh2/h**2, och2/h**2, omegal
         !  write (*,*) 100._dl*h, temp0k, ypxe
         !end if

         call sdratio_length_calc(thryflag,difflenc,ratioarray(neffnum-nind+1,:))

         !if (neffval.eq.3._dl) then
         !  d_a = ang_diam(adcpl,hubparm)
         !  write (*,*) 'Angular Diameter:', d_a
         !end if


         !soundhor = shcalc(adcpl,neffval)
         soundhor = shcalc(adcod,neffval)

         ratiosd = soundhor/difflenc
         ratioarray(neffnum-nind+1,size(ratioarray,2)) = ratiosd
         neffarray(neffnum-nind+1) = neffval

         if (neffval.eq.3._dl) then
           ratiotest = ratiosd
           soundhortest = soundhor
           difflenctest = difflenc
         end if

#ifdef prllel
       if (rank.eq.0) then
#endif
         !write (*,*) neffval, ratiosd
         !write (104,*) neffval, ratiosd
         !write (*,*) neffval, ratiosd, soundhor, difflenc
#ifdef prllel
       end if !(rank.eq.0
#endif
       
100    end do


!------calculate numerical derivatives at endpoints.

       mirrorwidth = (neffrangeh - neffrangel)/real(neffnum, kind=dl)/10._dl

       neffend(1) = neffrangel - mirrorwidth
       if (neffend(1).le.0._dl) neffend(1) = 0._dl
       neffend(2) = neffrangel + mirrorwidth
       neffend(3) = neffrangeh - mirrorwidth
       neffend(4) = neffrangeh + mirrorwidth

       do nind=1,4 !end do at 200

         neffval = neffend(nind)

         !omeganu = omegag*7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval
         omeganu = hubparm%og*7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*neffval
         hubparm%on = omeganu
         !omegal = 1._dl - hubparm%om - omegag - omeganu !contribution from dark energy
         omegal = 1._dl - hubparm%om - hubparm%og - omeganu !contribution from dark energy
         !energy density of dark energy at current epoch in units of MeV**4:
         rhode = omegal*hubparm%rhoc0
         hubparm%ol = omegal

         call sdratio_length_calc(thryflag,difflenc,ratioend(nind,:))       

         soundhor = shcalc(adcod,neffval)

         !ratioend(nind) = soundhor/difflenc
         ratioend(nind,size(ratioend,2)) = soundhor/difflenc
200    end do

       dratio(:) = ratioend(2,:) - ratioend(1,:)
       do i=1,size(dratio, 1)
         if (dratio(i).ne.0._dl) then
           dneff2(i) = (neffend(2) - neffend(1))/dratio(i)
         else
           dneff2(i) = 0._dl
         end if
       end do

       dratio(:) = ratioend(4,:) - ratioend(3,:)
       do i=1,size(dratio, 1)
         if (dratio(i).ne.0._dl) then
           dneff1(i) = (neffend(4) - neffend(3))/dratio(i)
         else
           dneff1(i) = 0._dl
         end if
       end do



!------interpolate to find neff value.


       ratiosd = ratiosdcalc

       !call interpneff(neffnum,ratioarray,neffarray,ratiosd &
       !                  ,tneffpred,dneff1,dneff2)
       do i=1,size(ratiocalc)
         call interpneff(ratioarray(:,i),neffarray,ratiocalc(i) &
                         ,tneffpred,dneff1(i),dneff2(i))
         !call cubic_spline(ratioarray(:,i),neffarray,ratioddot(:,i),dneff1(i),dneff2(i))
         !call cubic_splint(ratioarray(:,i),neffarray,ratioddot(:,i),ratiocalc(i),tneffpred)
         if (neffevflag) write (96,*) real(i, kind=dl)/real((intnum-1)/4, kind=dl)*adcod,ratiocalc(i),tneffpred
         if (i.eq.size(ratiocalc)) tneff = tneffpred
       end do

       !call interpneff(ratioarray(:,size(ratiocalc)),neffarray,1.04136_dl/0.161375_dl &
       !                  ,tneff,dneff1(size(ratiocalc)),dneff2(size(ratiocalc)))
       !write (*,*) 1.04136/0.161375, tneff

       !runarray(summnuind) = tneffpred

       tneff = tneffpred
       !write (90,*) xi(1), tneffpred
       !write (90,'(/)')
       if (neffevflag) then
         write (96,*) summnu, tneffpred
         write (96,'(/)')
       end if

       !write (107, *) obh2, ypxe, summnu, ratiosd, tneff 

#ifdef prllel
       if (rank.eq.0) then
#endif
       !write (*,*) 'recom is done'
       !write (*,*) obh2, ypri, ratiosd, tneff
       !write (*,*) tneff, ratiosd, soundhor, difflenc
       !write (*,*) omeganu*h**2
#ifdef prllel
       end if !(rank.eq.0)
#endif

       end do !summnuind



       deallocate(ratioarray)
       !deallocate(ratioddot)
       deallocate(ratioend)
       deallocate(ratiocalc)
       deallocate(dratio)
       deallocate(dneff1)
       deallocate(dneff2)
       deallocate(neffarray)
       deallocate(xec)
       !deallocate(xei)
       nullify(xei)
       deallocate(ascarray)
       !deallocate(tbar)
       nullify(tbar)
       !nullify(ratioy)
       deallocate(ratioy)
       deallocate(equilfrac)
       nullify(hubparm%nuparm)

 

       return


       end subroutine sdratio

