
       module xe_history


!------linkages.

!      used by - [subroutine] integ_a


!-----remarks.

!     Contains code to calculate the free-electron fraction using only hydrogen evolution.


!------modules.

       use sdratiovar


       implicit none


!------conversion factors.

       real(dl) :: convfact3 = 1._dl/2.99792458e+10_dl/(197.3269718*1.e-13_dl)**2
       real(dl) :: convfact4 = 1._dl/2.99792458e+08_dl/(197.3269718*1.e-15_dl)**2


!------transition scale factors.

       real(dl) atranshe32 !transition between HeIII and HeII
       real(dl) atranshe21 !transition between HeII and HeI
       real(dl) atransh21 !transition between HII and HI

       real(dl) neffvalxe
       integer lastind

       real(dl) bor_xe_de_h, bor_xe_de_he

       contains


!--------------------------------------------------------------------

       subroutine get_sel_recom_history(ratioy,yi,ascarray,astart,aend,fixflag)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Computes recombination history over a specified epoch.


       implicit none


       save


!------throughput variables.

       real(dl), dimension(:,:), intent(out) :: ratioy !array of evolution variables
       real(dl), dimension(:), intent(in) :: yi !initial values of evolution variables
       real(dl), dimension(:), intent(out) :: ascarray !array of scale-factor values
       real(dl), intent(in) :: astart
       real(dl), intent(in) :: aend
       logical, intent(in) :: fixflag !flag to use fixed stepsize


!------local variables.

       real(dl) apri, asc, astep, dasc !scale-factor variables
       integer i, j, intnum2


!------procedure.

       !write (*,*) hubparm%rhoc0, hubparm%og, hubparm%om, hubparm%ol, hubparm%orad, hubparm%temp0ev, hubparm%on, hubparm%onflag
       !write (*,*) hubparm%nuparm(1,1)%occfrac(1:3)
       !write (*,*) hubparm%nuparm(:,:)%mass
       !write (*,*) hubparm%nuparm(:,:)%dilfact
       !dasc = aend/real(intnum-1, kind=dl)
       dasc = (aend - astart)/real(intnum-1, kind=dl)

       intnum2 = intnum
       !if ((neffval.eq.3._dl).and.(xeflag)) intnum2 = pastfact*intnum
       if ((.not.fixflag).and.(xeflag)) intnum2 = pastfact*intnum
       !dasc = (aend - astart)/real(intnum-1, kind=dl)
       dasc = (aend - astart)/real(intnum2-1, kind=dl)
       !dascdy = dasc
       if (fixflag) then
         astep = dasc
       else
         astep = dasc/100._dl
       end if
       asc = astart - astep

       !do i=1,intnum
       do i=1,intnum2
         if (fixflag) then
           !asc = real(i-1, kind=dl)/real(intnum-1, kind=dl)*(aend - astart) + astart
           asc = real(i-1, kind=dl)/real(intnum2-1, kind=dl)*(aend - astart) + astart
         else
           astep = min(astep, (aend - asc))
           asc = asc + astep
         end if
         ascarray(i) = asc
         if (i.eq.1) then
           ratioy(i,1) = yi(1)
           ratioy(i,2) = yi(2)
           ratioy(i,3) = yi(3)
           ratioy(i,4) = yi(4)
         else
           !apri = real(i-2, kind=dl)/real(intnum-1, kind=dl)*(aend - astart) + astart
           !apri = real(i-2, kind=dl)/real(intnum2-1, kind=dl)*(aend - astart) + astart
           apri = asc - astep
           ratioy(i,:) = ratioy(i-1,:) !set new value of xei to old, and then proceed to integrate if applicable
           !astep = dasc
           if (asc.lt.atranshe21) then
             call saha_net(asc,ratioy(i,1:3))
             ratioy(i,4) = t0mev/asc
             !write (*,*) 'here3', ratioy(i,4)
           !else if (asc.lt.atransh21) then
           else if (ratioy(i,1).gt.(0.99_dl*xeimax(1))) then
             !Boltzmann integration for singly-ionized He:
             !call xe_he_driver(apri,3,xei(i,:),astep)
             call xe_driver(apri,ratioy(i,:),astep,recom_he_net,fixflag,.true.)
             !Input Saha relation for singly-ionized H
             !xei(i,1) = saha_ion((1._dl - ypxe/2._dl),asc,1._dl,ipoth1,xeimax(1))
           else
             !Boltzmann integration:
             !call xe_full_driver(apri,3,xei(i,:),astep)
             call xe_driver(apri,ratioy(i,:),astep,recom_full_net,fixflag,.false.)
           end if
         end if
         !if (.not.fixflag) then
         !  write (101,*) apri, asc, astep, sum(xei(i,:)*xeiz(:)), (xei(i,j),j=1,3)
         !end if
         !if (.not.fixflag) write (103,*) asc, i
         if (asc.ge.aend) then
           lastind = i
           !if (.not.fixflag) write (*,*) 'ended', lastind, asc, aend, astep
           exit
         end if
       end do


       return


       end subroutine get_sel_recom_history


!--------------------------------------------------------------------

       subroutine saha_net(asc,xei)


!------linkages.

!      called by - [subroutine] recom


!------remarks.

!      Computes recombination history over a specified epoch.


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: asc !scale factor
       real(dl), dimension(:), intent(inout) :: xei !array of free-electron fractions


!------procedure.

       if (asc.lt.atranshe32) then
         !Input Saha relation for doubly-ionized He:
         !xei(3) = saha_ion(ypxe/2._dl,asc,1._dl,ipothe2,xeimax(3))
         !xei(3) = saha_ion(ypxe/4._dl,asc,1._dl,ipothe2,xeimax(3))
         xei(3) = saha_ion(asc,1._dl,ipothe2,xeimax(3),xeimax(3))
         if (xei(3).le.(1.e-04_dl*xeimax(3))) then
           xei(3) = 0._dl
           xei(2) = xeimax(2) !No neutral He
         else
           xei(2) =  xeimax(2) - xei(3) !No neutral He
         end if
         xei(1) = xeimax(1)
       else if (asc.lt.atranshe21) then
         xei(3) = 0._dl
         !Input Saha relation for singly-ionized He:
         !xei(2) = saha_ion(ypxe/2._dl,asc,4._dl,ipothe1,xeimax(2))
         xei(2) = saha_ion(asc,4._dl,ipothe1,2._dl*xeimax(2),xeimax(2))
         xei(1) = xeimax(1)
       else !should not use this else statement - use Boltzmann equation integration
         xei(3) = 0._dl
         !xei(2) = saha_ion(ypxe/2._dl,asc,4._dl,ipothe1,xeimax(2))
         xei(2) = saha_ion(asc,4._dl,ipothe1,(2._dl*xeimax(2) + xei(1)),xeimax(2))
         !xei(1) = saha_ion((1._dl - ypxe/2._dl),asc,1._dl,ipoth1,xeimax(1))
         xei(1) = saha_ion(asc,1._dl,ipoth1,(xeimax(1) + 2._dl*xeimax(2) - xei(2)),xeimax(1))
       end if


       return


       end subroutine saha_net


!--------------------------------------------------------------------

       !subroutine xe_driver(apri,n,xe,astep0,f,fixflag,hflag)
       subroutine xe_driver(apri,ratioy,astep0,f,fixflag,hflag)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Computes free-electron fraction for a network.


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: apri !scale factor
       real(dl), dimension(:), intent(inout) :: ratioy !array of evolution variables
       real(dl), intent(inout) :: astep0 !step size in scale factor
       !real(dl), intent(in) :: tol !error tolerance
       logical, intent(in) :: fixflag !flag to not use dynamical timestep
       logical, intent(in) :: hflag !flag to calculate H by Saha
       external f


!------local variables.

       integer i, finegrid
       real(dl), dimension(:), allocatable :: yscal !array of error values for y
       real(dl) astep !current step size in scale factor
       real(dl) asteptemp !current step size in scale factor
       real(dl) anext !target value of scale factor needed to exit subroutine
       real(dl), dimension(:), allocatable :: ytemp !temp. array of evolution variables at scale factor
       real(dl), dimension(:), allocatable :: yerr !temp. array of free-electron fractions at scale factor
       real(dl) errmax, a
       real(dl) tol
       real(dl) dasum
       !real(dl), dimension(3) :: dydx


!------procedure.

       allocate(yscal(size(ratioy)))
       allocate(ytemp(size(ratioy)))
       allocate(yerr(size(ratioy)))

       tol = 1.e-07_dl
       finegrid = 100
       !astep = astep0
       a = apri
       if (fixflag) then
         astep = astep0/real(finegrid, kind=dl)
       else
         astep = astep0
       end if
       anext = a + astep0

       !yscal(:) = xeimax(:)
       yscal(1:3) = xeimax(1:3)
       yscal(4) = ratioy(4)
       ytemp = ratioy
       dasum = 0._dl
       do
         yerr = 0._dl
         call rk6_stepper(a,ytemp,f,astep,yerr)
         errmax = 0._dl
         do i=1,size(ratioy)
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
         end do
         errmax = errmax/tol
         if (errmax.gt.1._dl) then !did not meet error tolerance; redo
           ytemp = ratioy
           asteptemp = 0.9_dl*astep/errmax**0.25_dl
           astep = max(asteptemp, 0.1_dl*astep)
           astep = min(astep, (anext - a))
           if (astep.lt.(1.e-08_dl*astep0)) then
             write (*,*) 'issue in xe_driver'
             if (neffvalxe.eq.5._dl) write (100, *) 'issue in xe_driver'
             exit
           end if
         else !met error tolerance
           dasum = dasum + astep
           a = a + astep
           ratioy = ytemp
           if (errmax.gt.1.89e-04_dl) then
             astep = 0.9_dl*astep/errmax**0.2_dl
           else
             astep = 5._dl*astep
           end if
           if ((abs(dasum - astep0)/astep0).lt.1.e-04_dl) exit
           astep = min(astep, (anext - a))
         end if
       end do

       if (hflag) then
         !ratioy(1) = saha_ion((1._dl - ypxe/2._dl),a,1._dl,ipoth1,xeimax(1))
         !ratioy(1) = saha_ion((1._dl - ypxe),a,1._dl,ipoth1,xeimax(1))
         ratioy(1) = saha_ion(a,1._dl,ipoth1,(xeimax(1) + 2._dl*xeimax(2) - 0.016_dl),xeimax(1))
       end if

       if (.not.fixflag) astep0 = astep

       deallocate(yscal)
       deallocate(ytemp)
       deallocate(yerr)


       return


       end subroutine xe_driver


!--------------------------------------------------------------------

       subroutine xe_full_driver(apri,n,xe,astep0)
       !subroutine xe_full_driver(a,n,xe,astep,tol)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Computes free-electron fraction for full network.


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: apri !scale factor
       integer, intent(in) :: n !Number of ions
       real(dl), dimension(n), intent(inout) :: xe !array of free-electron fractions at scale factor
       !real(dl), dimension(:), intent(inout) :: xe !array of free-electron fractions at scale factor
       real(dl), intent(in) :: astep0 !step size in scale factor
       !real(dl), intent(in) :: tol !error tolerance


!------local variables.

       integer i, finegrid
       real(dl), dimension(n) :: yscal !array of error values for y
       real(dl) astep !current step size in scale factor
       real(dl) asteptemp !current step size in scale factor
       real(dl) anext !target value of scale factor needed to exit subroutine
       real(dl), dimension(n) :: xetemp !temp. array of free-electron fractions at scale factor
       real(dl), dimension(n) :: yerr !temp. array of free-electron fractions at scale factor
       real(dl) errmax, a
       real(dl) tol
       real(dl) dasum


!------procedure.

       tol = 1.e-07_dl
       finegrid = 100
       !astep = astep0
       a = apri
       astep = astep0/real(finegrid, kind=dl)
       anext = a + astep0

       !call rk4_stepper(a,n,xe,recom_full_net,astep,finegrid)
       !call rk6_stepper(a,n,xe,recom_full_net,astep,finegrid)
       yscal(:) = xeimax(:)
       xetemp = xe
       dasum = 0._dl
       do
         yerr = 0._dl
         !call rk6_stepper(a,n,xetemp,recom_full_net,astep,yerr)
         call rk6_stepper(a,xetemp,recom_full_net,astep,yerr)
         errmax = 0._dl
         do i=1,size(xe)
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
         end do
         errmax = errmax/tol
         if (errmax.gt.1._dl) then !did not meet error tolerance; redo
           xetemp = xe
           asteptemp = 0.9_dl*astep/errmax**0.25_dl
           astep = max(asteptemp, 0.1_dl*astep)
           astep = min(astep, (anext - a))
           if (astep.lt.(1.e-08_dl*astep0)) then
             write (*,*) 'issue in xe_full_driver'
             exit
           end if
         else !met error tolerance
           dasum = dasum + astep
           !if (neffvalxe.eq.3._dl) write (100, *) a, astep, astep0, dasum, anext, 'success'
           a = a + astep
           xe = xetemp
           if ((abs(dasum - astep0)/astep0).lt.1.e-04_dl) exit
           if (errmax.gt.1.89e-04_dl) then
             astep = 0.9_dl*astep/errmax**0.2_dl
           else
             astep = 5._dl*astep
           end if
           astep = min(astep, (anext - a))
         end if
       end do


       return


       end subroutine xe_full_driver


!--------------------------------------------------------------------

       subroutine xe_he_driver(a,n,xe,astep)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Computes free-electron fraction only for He.


       implicit none


       save


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       integer, intent(in) :: n !Number of ions
       real(dl), dimension(n), intent(inout) :: xe !array of free-electron fractions at scale factor
       !real(dl), dimension(:), intent(inout) :: xe !array of free-electron fractions at scale factor
       real(dl), intent(in) :: astep !step size in scale factor


!------local variables.

       integer finegrid


!------procedure.

       finegrid = 100

       call rk4_stepper(a,n,xe,recom_he_net,astep,finegrid)
       !call rk6_stepper(a,n,xe,recom_he_net,astep,finegrid)


       return


       end subroutine xe_he_driver


!--------------------------------------------------------------------

       subroutine recom_full_net(a,yin,dyda)


!------linkages.

!      called by - [function] rk6_stepper, get_xe
!      calls     - [function] hubcalc, alpha2, beta


!------remarks.

!      Subroutine to control network free-electron fraction X_e and baryon temperature tbar.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), dimension(:), intent(in) :: yin !evolution variables
       real(dl), dimension(:), intent(out) :: dyda !derivatives of yin


!------local variables.


!------procedure.

       !write (*,*) 'here1', size(yin)
       !HII:
       dyda(1) = xe_de_h(a,yin)
       bor_xe_de_h = dyda(1)
       !HeII:
       dyda(2) = xe_de_he(a,yin)
       bor_xe_de_he = dyda(2)
       !HeIII:
       dyda(3) = 0._dl
       !tbar:
       !dyda(4) = 0._dl
       dyda(4) = tbar_de(a,yin)


       return


       end subroutine recom_full_net


!--------------------------------------------------------------------

       subroutine recom_he_net(a,yin,dyda)


!------linkages.

!      called by - [function] rk6_stepper, get_xe
!      calls     - [function] hubcalc, alpha2, beta


!------remarks.

!      Subroutine to control He free-electron fraction X_e.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), dimension(:), intent(in) :: yin !evoltion variables
       real(dl), dimension(:), intent(out) :: dyda !derivatives of yin


!------local variables.


!------procedure.

       !HII:
       dyda(1) = 0._dl
       bor_xe_de_h = 0._dl
       !HeII:
       dyda(2) = xe_de_he(a,yin)
       bor_xe_de_h = dyda(2)
       !HeIII:
       dyda(3) = 0._dl
       !tbar:
       !dyda(4) = 0._dl
       dyda(4) = tbar_de(a,yin)


       return


       end subroutine recom_he_net


!--------------------------------------------------------------------

       real(dl) function saha_ion(a,spinfac,ipot,c1,c2)
       

!------linkages.

!      called by - [subroutine] get_xe_history, saha_net
!      calls     - none


!------remarks.

!      Calculates ionization fraction assuming Saha equilibrium.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), intent(in) :: spinfac !ratio of total spin quantum numbers
       real(dl), intent(in) :: ipot !ionization potential of ion (MeV)
       real(dl), intent(in) :: c1 !First constant needed
       real(dl), intent(in) :: c2 !Second constant needed


!------local variables.

       real(dl) tempmev !Temp. in MeV
       real(dl) ndens !Number density in MeV^3
       real(dl) rhs !Right-hand-side of Saha Eqn.
       real(dl) bcoeff !b coefficient of quadratic eqn.


!------procedure.

       !Temperature in MeV:
       tempmev = temp0ev*1.e-06_dl*(a0/a)

       !Number density of electrons at scale factor a:
       ndens = (1._dl - ypxe/2._dl)*obh2*3._dl/8._dl/pi/amu*mpls**2*convfact2*(a0/a)**3

       !Right-hand-side of Saha eqn.:
       rhs = spinfac/ndens*exp(-ipot/tempmev)*(melec*tempmev/2._dl/pi)**1.5_dl

       bcoeff = rhs + 1._dl - c1
       if (bcoeff.gt.1.e+08_dl) then !due to inaccuracy in double precision
         saha_ion = (1._dl - c2/rhs)*c2 !small value approximation
       else
         saha_ion = bcoeff/2._dl*(sqrt(1._dl + 4._dl*c2*rhs/bcoeff**2) - 1._dl) !quadratic formula
       end if
      

       return


       end function saha_ion


!--------------------------------------------------------------------

       !real(dl) function dsida(abund,a,spinfac,ipot,abundmax)
       real(dl) function dsida(a,spinfac,ipot,c1,c2)
       

!------linkages.

!      called by - [subroutine] get_xe_history
!      calls     - none


!------remarks.

!      Calculates derivative of saha_ion.


       implicit none


!------throughput variables.

       !real(dl), intent(in) :: abund !abundance of given nucleus
       real(dl), intent(in) :: a !scale factor
       real(dl), intent(in) :: spinfac !ratio of total spin quantum numbers
       real(dl), intent(in) :: ipot !ionization potential of ion (MeV)
       !real(dl), intent(in) :: abundmax !Maximum free-electrons for a given ion
       real(dl), intent(in) :: c1 !First constant needed
       real(dl), intent(in) :: c2 !Second constant needed


!------local variables.

       real(dl) tempmev !Temp. in MeV
       real(dl) ndens !Number density in MeV^3
       real(dl) rhs !Right-hand-side of Saha Eqn.
       real(dl) drhsda !Deriv. of right-hand-side of Saha Eqn.
       real(dl) bcoeff !term in quadratic formula
       real(dl) sqterm !square-root term


!------procedure.

       !Temperature in MeV:
       tempmev = temp0ev*1.e-06_dl*(a0/a)

       !Number density of ion of interest at scale factor a:
       !ndens = abund*obh2*3._dl/8._dl/pi/amu*mpls**2*convfact2*(a0/a)**3
       ndens = (1._dl - ypxe/2._dl)*obh2*3._dl/8._dl/pi/amu*mpls**2*convfact2*(a0/a)**3

       !Right-hand-side of Saha eqn.:
       rhs = spinfac/ndens*exp(-ipot/tempmev)*(melec*tempmev/2._dl/pi)**1.5_dl

       !Derivative of right-hand-side of Saha eqn.:
       drhsda = rhs/a*(1.5_dl - ipot/tempmev)

       !if (rhs.gt.1.e+08_dl) then !due to inaccuracy in double precision
       !  dsida = 1._dl/rhs**2*drhsda*abundmax !from small value approximation
       !else
       !  dsida = drhsda/sqrt(1._dl + 4._dl/rhs)/2._dl*abundmax* &
       !          (1._dl + 2._dl/rhs - sqrt(1._dl + 4._dl/rhs)) !from quadratic formula
       !end if

       bcoeff = rhs + 1._dl - c1
       if (bcoeff.gt.1.e+08_dl) then !due to inaccuracy in double precision
         dsida = c2**2/rhs**2*drhsda !small value approximation
       else
         sqterm = sqrt(1._dl + 4._dl*c2*rhs/bcoeff**2)
         dsida = drhsda/sqterm/2._dl*(1._dl + 2._dl*c2/bcoeff - sqterm) !from quadratic formula
       end if

      
       return


       end function dsida


!-----------------------------------------------------------

       !real(dl) function sinrsolver(f,dfdx,sol,tol,x,abund,spinfac,ipot,abundmax) !Newton-Raphson solver
       real(dl) function sinrsolver(f,dfdx,sol,tol,x,spinfac,ipot,c1,c2) !Newton-Raphson solver

       real(dl), external :: f
       real(dl), external :: dfdx
       real(dl), intent(in) :: sol !solution
       real(dl), intent(in) :: tol !tolerance
       real(dl), intent(in) :: x !input guess
       !real(dl), intent(in) :: abund !input into f and dfdx
       real(dl), intent(in) :: spinfac !input into f and dfdx
       real(dl), intent(in) :: ipot !input into f and dfdx
       !real(dl), intent(in) :: abundmax !input into f and dfdx
       real(dl), intent(in) :: c1 !First constant needed
       real(dl), intent(in) :: c2 !Second constant needed

       integer counter, counter2, counter3
       real(dl) xit

       counter = 1
       counter2 = 0
       counter3 = 0
       xit = x
       do
         !if (abs(f(abund,xit,spinfac,ipot,abundmax) - sol).lt.tol) exit
         if (abs(f(xit,spinfac,ipot,c1,c2) - sol).lt.tol) exit
         if (counter.eq.50) exit
         if (counter2.eq.5) exit
         !if (dfdx(abund,xit,spinfac,ipot,abundmax).ne.0._dl) then
         if (dfdx(xit,spinfac,ipot,c1,c2).ne.0._dl) then
           !write (*,*) xit, f(abund,xit,spinfac,ipot,abundmax), dfdx(abund,xit,spinfac,ipot,abundmax)
           !xit = xit - (f(abund,xit,spinfac,ipot,abundmax) - sol)/dfdx(abund,xit,spinfac,ipot,abundmax)
           xit = xit - (f(xit,spinfac,ipot,c1,c2) - sol)/dfdx(xit,spinfac,ipot,c1,c2)
           if (xit.le.0._dl) then
             xit = x + 1.e-05_dl*real(counter3 + 1, kind=dl)
             counter3 = counter3 + 1
           end if
         else
           xit = x + 1.e-05_dl*real(counter2 + 1, kind=dl)
           counter2 = counter2 + 1
           counter = 0
           write (*,*) 'Bad input into nrsolver'
         end if
         counter = counter + 1
       end do

       sinrsolver = xit


       end function sinrsolver


!--------------------------------------------------------------------

       real(dl) function xe_de_h(a,yin)
       

!------linkages.

!      called by - [subroutine] recom_full_net
!      calls     - [function] hubcalc, alphab, beta


!------remarks.

!      Differential equation to compute free-electron fraction X_e^{(H)}.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), dimension(:), intent(in) :: yin !evolution variables


!------local variables.

       real(dl) x
       real(dl) xeh
       real(dl) xetot
       real(dl) tbar !baryon temp.
       real(dl) corr_fact
       real(dl) lama
       real(dl) beta2
       real(dl) :: lamgam = 8.2245809_dl*6.58211928e-22_dl !Rate for 2-photon emission, units of MeV
       real(dl) expfact !factor to avoid overflow error


!------procedure.

       xetot = sum(yin(1:3)*xeiz(1:3)) !Total free electron fraction
       xeh = yin(1) !Hydrogen free-electron fraction
       tbar = yin(4) !baryon temp.
       !x = ipoth1/t0mev*a !Integration variable x
       x = ipoth1/tbar !Integration variable x
       
       if (x.eq.0._dl) then
         xe_de_h = 0._dl
       else
         if (xeh.ne.xeimax(1)) then
           if ((0.75*x).lt.300._dl) then
             expfact = exp(0.75_dl*x)
           else
             expfact = exp(300._dl)
           end if
           !beta2 = alphab(t0mev/a)*beta(t0mev/a,ipoth1)*expfact !Ionization coefficient
           beta2 = alphab(tbar)*beta(tbar,ipoth1)*expfact !Ionization coefficient
           lama = hubcalc(a,hubparm)*(3._dl*ipoth1)**3/(8._dl*pi)**2 &
                  /3._dl*8._dl*pi/mpls**2/obh2*amu/(xeimax(1) - xeh) &
                  *a**3/convfact2 !Rate for Lyman-alpha emission
           corr_fact = (lama + lamgam)/(lama + lamgam + beta2) !Peebles correction factor
         else
           corr_fact = 1._dl
         end if
         !xe_de_h = alphab(t0mev/a)/hubcalc(a,hubparm)/a &
         xe_de_h = alphab(tbar)/hubcalc(a,hubparm)/a &
                   !*((xeimax(1) - xeh)*beta(t0mev/a,ipoth1) &
                   *((xeimax(1) - xeh)*beta(tbar,ipoth1) &
                   - xeh*xetot*(1._dl - ypxe/2._dl) &
                   *3._dl/8._dl/pi*mpls**2*obh2/amu*convfact2/a**3)*corr_fact
       end if
       
       !if (abs(xe_de_h*a).lt.1.e-05_dl) xe_de_h = 0._dl
       if (xe_de_h.gt.0._dl) xe_de_h = 0._dl


       return


       end function xe_de_h


!--------------------------------------------------------------------

      real(dl) function xe_de_he(a,yin)
       

!------linkages.

!      called by - [subroutine] recom_full_net, recom_he_net
!      calls     - [function] hubcalc, alpha2, beta


!------remarks.

!      Differential equation to compute free-electron fraction X_e^{(He)}.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), dimension(:), intent(in) :: yin !evolution variables


!------local variables.

       real(dl) xehe
       real(dl) xetot
       real(dl) tbar !baryon temperature
       real(dl) he1h !factors involving ypxe
       real(dl) corr_fact
       real(dl) lama
       real(dl) beta2
       real(dl) :: lamgam = 8.2245809_dl*6.58211928e-22_dl !Rate for 2-photon emission, units of MeV
       real(dl) :: lamgam2 = 51.3_dl*6.58211928e-22_dl !Rate for 2-photon emission, units of MeV
       real(dl) he_boltz !Factor for Peebles correction
       real(dl) expfact !factor to avoid overflow error


!------procedure.

       xetot = sum(yin(1:3)*xeiz(1:3)) !Total free electron fraction
       xehe = yin(2) !singlely ionized Helium free-electron fraction
       tbar = yin(4) !baryon temperature

       he1h = (1._dl - ypxe)/(1._dl - ypxe/2._dl)
       !6.02247886e-07 MeV = difference in energies between 2p and 2s levels in HeII

       if (a.eq.0._dl) then
         xe_de_he = 0._dl
       else
         if (xehe.ne.xeimax(2)) then
           !if ((he1c2s*a/t0mev).lt.300._dl) then
           !  expfact = exp(he1c2s*a/t0mev)
           if ((he1c2s/tbar).lt.300._dl) then
             expfact = exp(he1c2s/tbar)
           else
             expfact = exp(300._dl)
           end if
           !if ((6.02247886e-07_dl*a/t0mev).lt.680._dl) then
           !  he_boltz = exp(6.02247886e-07_dl*a/t0mev)
           if ((6.02247886e-07_dl/tbar).lt.680._dl) then
             he_boltz = exp(6.02247886e-07_dl/tbar)
           else
             he_boltz = exp(680._dl)
           end if
           !beta2 = alphahe(t0mev/a)*4._dl*beta(t0mev/a,ipothe1)*expfact !Ionization coefficient
           beta2 = alphahe(tbar)*4._dl*beta(tbar,ipothe1)*expfact !Ionization coefficient
           lama = hubcalc(a,hubparm)*(0.86296316268_dl*ipothe1)**3/pi**2 &
                  /3._dl*8._dl*pi/mpls**2/obh2*amu/(xeimax(2) - xehe) &
                  *a**3/convfact2 !Rate for Lyman-alpha emission
           !Problem here: (02 Jan 2015)
           corr_fact = (lama + lamgam2*he_boltz)/(lama + (lamgam2 + beta2)*he_boltz) !Peebles correction factor
         else
           corr_fact = 1._dl
         end if
         !xe_de_he = alphahe(t0mev/a)/hubcalc(a,hubparm)/a &
         !          *((xeimax(2) - xehe)*4._dl*beta(t0mev/a,ipothe1) & !factor of 4 from stat. weights
         xe_de_he = alphahe(tbar)/hubcalc(a,hubparm)/a &
                   *((xeimax(2) - xehe)*4._dl*beta(tbar,ipothe1) & !factor of 4 from stat. weights
                   - xehe*xetot*(1._dl - ypxe/2._dl) &
                   *3._dl/8._dl/pi*mpls**2*obh2/amu*convfact2/a**3)*corr_fact
       end if
       
       !if (abs(xe_de_he*a).lt.1.e-05_dl) xe_de_he = 0._dl
       if (xe_de_he.gt.0._dl) xe_de_he = 0._dl

       
       !if (neffvalxe.eq.5._dl) write (100, *) corr_fact, lama, he_boltz, beta2, 6.02247886e-07_dl*a/t0mev

       return


       end function xe_de_he


!--------------------------------------------------------------------

       real(dl) function alpha2(x)
       

!------linkages.

!      called by - [function] xe_de_h
!      calls     - none


!------remarks.

!      Computes H recombination coefficient.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x


!------procedure.

       if (x.lt.1._dl) then
         alpha2 =  0._dl
       else
         alpha2 = 9.78_dl*alpha**2/melec**2*sqrt(x)*log(x)
       end if
       

       return


       end function alpha2


!--------------------------------------------------------------------

       real(dl) function alphab(tempmev)
       

!------linkages.

!      called by - [function] xe_de_h
!      calls     - none


!------remarks.

!      Computes H recombination coefficient.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: tempmev


!------local variables.

       real(dl) tk !Temperature in K


!------procedure.

       tk = tempmev/boltzmann*1.e+15_dl

       alphab = 1.e-13_dl*4.309_dl*(1.e-04_dl*tk)**(-0.6166_dl) &
                /(1._dl + 0.6703_dl*(1.e-04_dl*tk)**0.5300_dl) &
                *convfact3
       

       return


       end function alphab


!--------------------------------------------------------------------
       
       real(dl) function beta(tempmev,ipot)


!------linkages.

!      called by - [function] xe_de_c
!      calls     - none


!------remarks.

!      Used for H ionization coefficient.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: tempmev !temperature in MeV
       real(dl), intent(in) :: ipot !ionization potential


!------procedure.
       
       if (tempmev.gt.1.e+30_dl) then
         beta = 0._dl
       else
         beta = (melec*tempmev/2._dl/pi)**1.5_dl*exp(-ipot/tempmev)
       end if
       
       return


       end function beta


!--------------------------------------------------------------------

       real(dl) function alphahe(tempmev)
       

!------linkages.

!      called by - [function] xe_de_he
!      calls     - none


!------remarks.

!      Computes recombination coefficient onto neutral He.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: tempmev


!------local variables.

       real(dl) tk !Temperature in K
       real(dl) term1, term2, term3 !terms needed for coefficient


!------procedure.

       tk = tempmev/boltzmann*1.e+15_dl

       term1 = sqrt(tk/3._dl)
       term2 = (1._dl + sqrt(tk/3._dl))**(1._dl - 0.711_dl)
       term3 = (1._dl + sqrt(tk/10._dl**(5.114_dl)))**(1._dl + 0.711_dl)

       alphahe = 10._dl**(-16.744_dl)/(term1*term2*term3)*convfact4
       

       return


       end function alphahe


!--------------------------------------------------------------------

      real(dl) function tbar_de(a,yin)
       

!------linkages.

!      called by - [subroutine] recom_full_net, recom_he_net
!      calls     - none


!------remarks.

!      Differential equation to compute baryon temperature.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: a !scale factor
       real(dl), dimension(:), intent(in) :: yin !evolution variables


!------local variables.

       real(dl) tbar !baryon temperature
       real(dl) timeh, timeth !time scale estimates 
       real(dl) xec !ionization fraction
       real(dl) prefac !term in smoothing function
       real(dl) term1, term2, term3, term4 !terms involving number density ratios


!------procedure.

       tbar = yin(4)

       timeh = 2._dl/3._dl/hubcalc(a,hubparm)

       xec = sum(yin(1:3)*xeiz(1:3))

       term1 = xec*(1._dl - ypxe/2._dl) + 1._dl - 3._dl*ypxe/4._dl
       timeth = 3._dl*melec/8._dl/sigmat*(15._dl/pi**2*(a/t0mev)**4)

       if (timeth.lt.(1.e-03_dl*timeh)) then
         !use smoothing term (arXiv:0902.3438)
         term2 = term1/(1._dl - 3._dl*ypxe/4._dl)
         term3 = 1._dl - ypxe
         !term3 = 1._dl
         term4 = ypxe/4._dl
         !term4 = ypxe/4._dl/(1._dl - ypxe)
         prefac = hubcalc(a,hubparm)/term1*timeth*t0mev/a
         !tbar_de = -1._dl/a**2*(t0mev &
         tbar_de = -t0mev/a**2 - prefac*3._dl/a &
                   + prefac/term2*(term3*bor_xe_de_h + term4*bor_xe_de_he)/xec &
                   - prefac*dhubda(a,hubparm)/hubcalc(a,hubparm)
         !tbar_de = -tbar/a
       else
         !use derivative
         tbar_de = -1._dl/a*((tbar - t0mev/a)*xec/term1/timeth/hubcalc(a, hubparm) + &
                   2._dl*tbar)
         !tbar_de = -2._dl*tbar/a
       end if

       !if (neffvalxe.eq.5._dl) write (102,*) a, xec, tbar_de, tbar, t0mev/a
       !if (neffvalxe.eq.5._dl) write (102,*) a,term1, term2, term3, term4

       !tbar_de = -tbar/a

       return


       end function tbar_de


       end module xe_history
