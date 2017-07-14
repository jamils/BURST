
       subroutine sdratio_init
       
!------linkages.

!      called by - [subroutine] main_init
!      calls     - none


!------remarks.

!      Initiates sdratio variables.


!------modules.

       use sdratiovar
       !allocatable arrays live here:
       use sdratio_length


       implicit none


!------local variables.


       real(dl) omegag !contribution from photons
       real(dl) rhoc0


!------data file.


10     format (6x,'neff',8x,'kdc',9x,'rsi',9x,'rsx',8x,'ratio',/,60('-'))

25     format (1p10e12.3)



!------procedure.

!------set values from sdratio_params.ini.

       open (unit=30, file='sdratio_params.ini', status='unknown')

       read (30,*) reddcpl
       read (30,*) neffnum
       read (30,*) neffwidth
       read (30,*) intnum
       read (30,*) rk4num
       read (30,*) summnuprec
       read (30,*) summnulow
       read (30,*) summnuhigh
       read (30,*) hierflag
       read (30,*) masstol
       read (30,*) xeflag
       read (30,*) pastfact
       read (30,*) neffevflag

       close (unit=30)


       if (mod(intnum,4).ne.1) then
         intnum = 4*(intnum/4) + 1
       end if

       neffevnum = (intnum-1)/4

       if (neffevflag) then
         allocate(ratioarray(neffnum,neffevnum)) !assign dimension of ratioarray
         !allocate(ratioddot(neffnum,neffevnum)) !assign dimension of ratioddot
         allocate(ratioend(4,neffevnum)) !assign dimension of ratioend
         allocate(ratiocalc(neffevnum)) !assign dimension of ratiocalc
         allocate(dratio(neffevnum)) !assign dimension of dratio
         allocate(dneff1(neffevnum)) !assign dimension of dneff1
         allocate(dneff2(neffevnum)) !assign dimension of dneff2
       else
         allocate(ratioarray(neffnum,1)) !assign dimension of ratioarray
         !allocate(ratioddot(neffnum,1)) !assign dimension of ratioddot
         allocate(ratioend(4,1)) !assign dimension of ratioend
         allocate(ratiocalc(1)) !assign dimension of ratiocalc
         allocate(dratio(1)) !assign dimension of dratio
         allocate(dneff1(1)) !assign dimension of dneff1
         allocate(dneff2(1)) !assign dimension of dneff2
       end if
       allocate(neffarray(neffnum)) !assign dimension of neffarray
       if (xeflag) then
         allocate(xec(pastfact*intnum)) !assign dimension of xec
         !allocate(xei(pastfact*intnum,3)) !assign dimension of xei
         allocate(ascarray(pastfact*intnum)) !assign dimension of ascarray
         !allocate(tbar(pastfact*intnum)) !assign dimension of tbar
         allocate(ratioy(pastfact*intnum,4)) !assign dimension of ratioy
       else
         allocate(xec(intnum)) !assign dimension of xec
         !allocate(xei(intnum,3)) !assign dimension of xei
         allocate(ascarray(intnum)) !assign dimension of ascarray
         !allocate(tbar(intnum)) !assign dimension of tbar
         allocate(ratioy(intnum,4)) !assign dimension of ratioy
       end if
       allocate(equilfrac(nbins)) !assign dimension of equilfrac
       !allocate(equilfrac(64)) !assign dimension of equilfrac



       temp0ev = boltzmann*temp0k !photon temp. at current epoch in eV
       t0mev = temp0ev*1.e-06_dl !photon temp. at current epoch in MeV
       hubparm%temp0ev = temp0ev

       adcpl = 1._dl/(reddcpl + 1._dl) !scale factor at decoupling from Planck in Mpc
       a0 = 1._dl !scale factor at current epoch in Mpc

       !photon energy density at current epoch in MeV**4:
       rhog0 = 2._dl*pi**2/30._dl*t0mev**4
       !critical energy density at current epoch in MeV**4 (only place h used):
       rhoc0 = 3._dl/8._dl/pi*mpls**2*h**2*convfact2
       hubparm%rhoc0 = rhoc0

       omegag = rhog0/rhoc0 !contribution from photons
       hubparm%og = omegag

       !Contribution from dark radiation:
       hubparm%orad = 0._dl


 
       return


       end subroutine sdratio_init
