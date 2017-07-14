
       subroutine bbn_yes(obh2,neffparm,ypri,yd,yhe3,yli7,nubbn)
 
!------linkages.

!      called by - [program] main
!      calls     - [subroutine] bbn_no_step


!------remarks.

!      Computes transport without BBN.


!------modules.

       !use bbnvar
       use bbnvar_v2
       use msw


       implicit none


       interface bbn_yes_launch_v2_interface
         subroutine bbn_yes_launch_v2(t,dt,bbnvalues,obh2)
           use bbnvar_v2
           use bessel
           use nse
           use renorm
           use ratenp
           implicit none
           real(dl), intent(out) :: t !time
           real(dl), intent(out) :: dt !time step
           type(bbnevolvar), intent(out) :: bbnvalues !bbn evolution variable
           real(dl), intent(in) :: obh2 !Baryon number
         end subroutine bbn_yes_launch_v2
       end interface bbn_yes_launch_v2_interface

       interface bbn_yes_driver_v2_interface
         subroutine bbn_yes_driver_v2(t,bbnvalues,dt)
           use bbnvar_v2
           implicit none
           real(dl), intent(inout) :: t !time
           type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
           real(dl), intent(inout) :: dt !time step
         end subroutine bbn_yes_driver_v2
       end interface bbn_yes_driver_v2_interface

       interface bbn_yes_revert_v2_interface
         subroutine bbn_yes_revert_v2(t,dt,bbnvalues)
           use bbnvar_v2
           implicit none
           real(dl), intent(out) :: t !time
           real(dl), intent(out) :: dt !time step
           type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable
         end subroutine bbn_yes_revert_v2
       end interface bbn_yes_revert_v2_interface

       save


!------throughput variables.

       real(dl), intent(in) :: obh2 !Baryon number
       real(dl), intent(in) :: neffparm !Parameter for changing Neff
       real(dl), intent(out) :: ypri !4He Mass Fraction
       real(dl), intent(out) :: yd !D relative abundance
       real(dl), intent(out) :: yhe3 !3He relative abundance
       real(dl), intent(out) :: yli7 !7Li relative abundance
       type(nuvar), dimension(:,:), intent(inout), target :: nubbn !\nu variable


!------local variables.

       integer i,j,m,n !indicies
       integer inum                 !selection number.

       real(dl) t !time
       real(dl) dt !time step
       !type(bbnevolvar), pointer :: bbnvalues !bbn evolution variable
       type(bbnevolvar) :: bbnvalues !bbn evolution variable

       real(dl) sentar


!------procedure.


!------set values from throughput variables.

       eta1 = convfact3/(boltzmann*temp0k)**3*obh2
       sentar = convfact4*(boltzmann*temp0k*1.e-06_dl)**3/obh2

       if (neffparm.ge.0._dl) then
         rhors = neffparm
       else
         rhors = 0._dl
       end if


!------set values for individual trans run.

       allocate(bbnvalues%nuoccprob(size(nubbn,1),2)) !allocate size of bbnvalues%nuoccprob
       bbnvalues%nuoccprob => nubbn
       call bbn_yes_launch_v2(t,dt,bbnvalues,obh2) !input information for trans run.
       if (saveflag) call bbn_yes_revert_v2(t,dt,bbnvalues)
       !it = it + 1


!------launch MSW run

       if (mswflag) call msw_launch(bbnvalues%nuoccprob)


!------begin individual trans run.

       call bbn_yes_driver_v2(t,bbnvalues,dt) !do transport computation.
       

       do i = 1,it !Temperature in MeV.
         t9out(i) = t9out(i)*.08617
       end do


#ifdef prllel
       if (rank.eq.0) then
#endif

       if (bbnflag) then

       write (21,2000) runind, obh2
2000   format ('run: ',i6,'; obh2 = ',1pe9.3)

       write (21,2002) cy,ct,t9i,t9f,ytmin
2002   format (' computational parameters:',/, &
               '   cy = ',f5.3,'/  ct = ',f5.3, &
               '/  initial temp = ',1pe8.2, &
               '/  final temp = ',1pe8.2, &
               '/  smallest abundances allowed = ',1pe8.2)

       write (21,2004) g,mntau,3.0,0._dl,xi(1),xi(2),xi(3)
2004   format (' model parameters:',/, &
               '   g = ',f5.2,'/  tau = ',f6.2, &
               '/  # nu = ',f5.2,'/  lambda = ',1pe10.3, &
               '/  xi-e = ',e10.3,'/  xi-m = ',e10.3, &
               '/  xi-t = ',e10.3)

       write (21,2005) it
2005   format ('it = ',i6,/)

!..........print headings, abundances for neutron to li8.
       write (21,2006)
2006   format (4x,'temp',8x,'n/H',10x,'p',10x,'D/H',9x,'T/H',8x, &
                  '3He/H',8x,'4He',8x,'6Li/H',7x,'7Li/H',7x, &
                  '7Be/H',/,120('-'))
       do j = 1,it
         write (21,2008) t9out(j),(xout(j,i),i=1,9)
2008     format (1pe10.3,1p10e12.3)
       end do

!..........print thermodynamic quantities.
       write (21,2010)
2010   format (' ',/,4x,'temp',9x,'t',10x,'rhog',8x,'rhoe',7x, &
                        'rhone',8x,'rhob',8x,'phie',9x,'dt',9x, &
                        'sen',10x,'H',9x,'eta',/,132('-'))
       do j = 1,it
         write (21,2012) t9out(j),tout(j),(thmout(j,i),i=1,4), &
                         thmout(j,5),dtout(j),senout(j),hubout(j), &
                         thmout(j,8)
2012     format (1pe10.3,1p10e12.3)
       end do

       write (21,2014)
2014   format (///)

       end if !bbnflag

#ifdef prllel
       end if !(rank.eq.0)
#endif


       !ypri = y(6)*4._dl
       !yd = y(3)/y(2)
       !yhe3 = (y(4) + y(5))/y(2) !use A = 3 isobar
       !yli7 = (y(9) + y(8))/y(2) !use A = 7 isobar

       ypri = bbnvalues%y(6)*4._dl
       yd = bbnvalues%y(3)/bbnvalues%y(2)
       yhe3 = (bbnvalues%y(4) + bbnvalues%y(5))/bbnvalues%y(2) !use A = 3 isobar
       yli7 = (bbnvalues%y(9) + bbnvalues%y(8))/bbnvalues%y(2) !use A = 7 isobar

       !if (rank.eq.0) write(*,*) cy, is, ypri, yd, yhe3, yli7
       
       !deallocate(bbnvalues%nuoccprob)


       return


       end subroutine bbn_yes
