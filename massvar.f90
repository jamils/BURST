
       module massvar


!------remarks.

!      Module to control the lowest neutrino mass eigenstate.


!------modules.

       use mainvar


       implicit none


!------nu mixing angles (s = sin and c = cos).

       real(dl), parameter :: s23 = sqrt(0.386_dl)
       real(dl), parameter :: s13 = sqrt(0.0241_dl)
       real(dl), parameter :: s12 = sqrt(0.307_dl)

       real(dl), parameter :: c23 = sqrt(1._dl - s23**2)
       real(dl), parameter :: c13 = sqrt(1._dl - s13**2)
       real(dl), parameter :: c12 = sqrt(1._dl - s12**2)

!------nu mixing matrix (unitary transf. from weak to mass).

       real(dl), dimension(3,3) :: mmix = reshape([ &
                c12*c13 , s12*c23 - c12*s23*s13 , s12*s23 + c12*c23*s13 , &
               -s12*c13 , c12*c23 + s12*s23*s13 , c12*s23 - s12*c23*s13 , &
                 -s13   ,       -s23*c13        ,        c23*c13        ],[3,3])


       contains


!-----------------------------------------------------------

       subroutine get_numass(summnu,hierflag,numass,tol)


!------linkages.

!      called by - [program] main
!      calls     - [subroutine] mass_eigenstate

!------remarks.

!      Subroutine to calculate the neutrino mass eigenstates given a sum of the neutrino mass statistic.


       implicit none


!------throughput variables.

       real(dl), intent(inout) :: summnu !neutrino mass sum
       logical, intent(in) :: hierflag !Flag for normal (.true.) or inverted (.false.) hierarchy
       real(dl), dimension(3), intent(out) :: numass !nu mass values
       real(dl), intent(in) :: tol !Tolerance for accepting value in nrsolver


!------local variables.

       real(dl) massl !lowest neutrino mass eigenstate

!------procedure.

       if (summnu.eq.0._dl) then
         numass = 0._dl
       else
         call mass_eigenstate(summnu,hierflag,massl,tol)
         if (hierflag) then
           numass(1) = massl
           numass(2) = sqrt(numass(1)**2 + 7.54e-05_dl)
           numass(3) = sqrt(numass(1)**2 + 2.47e-03_dl) !numass(1) instead of numass(2)?
         else !inverted hiearchy
           numass(3) = massl
           numass(1) = sqrt(numass(3)**2 + 2.37e-03_dl)
           numass(2) = sqrt(numass(1)**2 + 7.54e-05_dl)
         end if
         summnu = sum(numass)
       end if


       return


       end subroutine get_numass


!-----------------------------------------------------------

       subroutine mass_eigenstate(summnu,hierflag,massl,tol)


!------linkages.

!      called by - [subroutine] recom
!      calls     - [function] nrsolver, fn, fi


!------remarks.

!      Subroutine to calculate the lowest neutrino mass eigenstate given a sum of the neutrino mass statistic.


       implicit none


       save



!------throughput variables.

       real(dl), intent(in) :: summnu !neutrino mass sum
       logical, intent(in) :: hierflag !Flag for normal (.true.) or inverted (.false.) hierarchy
       real(dl), intent(out) :: massl !Lowest mass eigenstate
       real(dl), intent(in) :: tol !Tolerance for accepting value in nrsolver


!-----procedure.

       if (hierflag) then
         if (fn(0._dl,summnu).ge.0._dl) then
           write (*,*) '\Sigma m_\nu too small,'
           write (*,*) 'setting smallest eigenstate to 1 meV'
           massl = 1.e-03_dl
         else
           massl = nrsolver(fn,dfndx,tol,summnu/3._dl,summnu)
         end if
       else
         if (fi(0._dl,summnu).ge.0._dl) then
           write (*,*) '\Sigma m_\nu too small,'
           write (*,*) 'setting smallest eigenstate to 1 meV,'
           massl = 1.e-03_dl
         else
           massl = nrsolver(fi,dfidx,tol,summnu/3._dl,summnu)
         end if
       end if


       return


       end subroutine mass_eigenstate


!-----------------------------------------------------------

       subroutine numass_assign(nuwm)


!------linkages.

!      called by - [program] main
!      calls     - none

!------remarks.

!      Subroutine to calculate the neutrino spectra in mass eigenstates.


       implicit none


!------throughput variables.

       type(nuvar), dimension(:,:), intent(inout) :: nuwm !\nu variable: weak basis -> mass basis


!------local variables.

       integer i,m,mm,n
       real(dl), dimension(:,:,:), allocatable :: tempoccfrac !nu occfrac in weak eigenbasis

!------procedure.

       allocate(tempoccfrac(3,2,size(nuwm(1,1)%occfrac))) !assign dimension of woccfrac

       do n=1,2
         do m=1,3
           tempoccfrac(m,n,:) = nuwm(m,n)%occfrac(:)
         end do !m
       end do !n
       do i=1,size(nuwm(1,1)%occfrac)
         do n=1,2
           do m=1,3 !for each mass eigenstate
             nuwm(m,n)%occfrac(i) = sum(mmix(m,:)**2*tempoccfrac(:,n,i))
           end do
         end do
       end do


       deallocate(tempoccfrac)


       return


       end subroutine numass_assign


!-----------------------------------------------------------

       real(dl) function nrsolver(f,dfdx,tol,x,summnu) !Newton-Raphson solver

       real(dl), external :: f
       real(dl), external :: dfdx
       real(dl), intent(in) :: tol !tolerance
       real(dl), intent(in) :: x !input guess
       real(dl), intent(in) :: summnu !sum of the neutrino mass

       integer counter
       real(dl) xit

       counter = 1
       xit = x
       do
         if (abs(f(xit,summnu)).lt.tol) exit
         if (counter.eq.50) exit
         if (dfdx(xit).ne.0._dl) then
           xit = xit - f(xit,summnu)/dfdx(xit)
         else
           xit = 1.e-03_dl
           write (*,*) 'Bad input into nrsolver'
           exit
         end if
         counter = counter + 1
       end do

       nrsolver = xit


       end function nrsolver


!-----------------------------------------------------------

       real(dl) function fn(x,summnu) !normal hierarchy

       real(dl), intent(in) :: x !input value
       real(dl), intent(in) :: summnu !input value for sum of the neutrino masses.


       fn = x + sqrt(x**2 + 7.54e-05_dl) + sqrt(x**2 + 2.47e-03_dl) - summnu


       end function fn


!-----------------------------------------------------------

       real(dl) function fi(x,summnu) !inverted hierarchy

       real(dl), intent(in) :: x !input value
       real(dl), intent(in) :: summnu !input value for sum of the neutrino masses.


       fi = x + sqrt(x**2 + 2.37e-03_dl) + sqrt(x**2 + 2.37e-03_dl + 7.54e-05_dl) - summnu


       end function fi


!-----------------------------------------------------------

       real(dl) function dfndx(x) !normal hierarchy derivative

       real(dl), intent(in) :: x !input value


       dfndx = 1._dl + x/sqrt(x**2 + 7.54e-05_dl) + x/sqrt(x**2 + 2.47e-03_dl)


       end function dfndx


!-----------------------------------------------------------

       real(dl) function dfidx(x) !inverted hierarchy derivative

       real(dl), intent(in) :: x !input value


       dfidx = 1._dl + x/sqrt(x**2 + 2.37e-03_dl) + x/sqrt(x**2 + 2.37e-03_dl + 7.54e-05_dl)


       end function dfidx


       end module massvar
