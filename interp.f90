
       module interp


!------linkages.

!      used by - [program] recom


!-----remarks.

!     Contains code to do interpolation for r_s/r_d


!------modules.

       use mainvar


       implicit none


       contains


!-----------------------------------------------------------

       subroutine interpneff(xvalues,yvalues,x,ypred,yp1,yp2)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Performs cubic spline and splint using points.


       implicit none


       save


!------throughput variables.

       real(dl), dimension(:), intent(in) :: xvalues
       real(dl), dimension(:), intent(in) :: yvalues
       real(dl), intent(in) :: x
       real(dl), intent(out) :: ypred
       real(dl), optional, intent(in) :: yp1
       real(dl), optional, intent(in) :: yp2


!------local variables.

       integer i,j
       integer eind

       integer arraysize !size of xvalues
       real(dl), dimension(:), allocatable :: yddot
       real(dl), dimension(:), allocatable :: ufact
       real(dl), dimension(:), allocatable :: vfact
       real(dl) wfact
       real(dl) sig1
       real(dl) sig2
       real(dl) para
       real(dl) parb
       real(dl) parh

!------procedure.


!------spline.

       arraysize = size(xvalues)

       allocate(yddot(arraysize)) !assign size of array yddot
       allocate(ufact(arraysize)) !assign size of array ufact
       allocate(vfact(arraysize)) !assign size of array vfact

       yddot = 0._dl
       ufact = 0._dl
       vfact = 0._dl
       wfact = 0._dl

       if (present(yp1)) then
         vfact(1) = -0.5_dl
         ufact(1) = 3._dl/(xvalues(2) - xvalues(1)) &
                    *((yvalues(2) - yvalues(1))/(xvalues(2) - xvalues(1)) - yp1)
       else
         vfact(1) = 0._dl
         ufact(1) = 0._dl
       end if
       
       do i=2,(arraysize-1)
         sig1 = (xvalues(i) - xvalues(i-1))/(xvalues(i+1) - xvalues(i-1))
         sig2 = sig1*vfact(i-1) + 2._dl
         vfact(i) = (sig1 - 1._dl)/sig2
         wfact = &
              (yvalues(i+1) - yvalues(i))/(xvalues(i+1) - xvalues(i)) &
              - (yvalues(i) - yvalues(i-1))/(xvalues(i) - xvalues(i-1))
         ufact(i) = (6._dl*wfact/(xvalues(i+1) - xvalues(i-1)) &
              - sig1*ufact(i-1))/sig2
       end do

       if (present(yp2)) then
         vfact(arraysize) = 0.5_dl
         ufact(arraysize) = 3._dl/(xvalues(arraysize) - xvalues(arraysize-1)) &
                            *(yp2 - (yvalues(arraysize) - yvalues(arraysize-1)) &
                            /(xvalues(arraysize) - xvalues(arraysize-1)))
       else
         vfact(arraysize) = 0._dl
         ufact(arraysize) = 0._dl
       end if
       
       yddot(arraysize) = (ufact(arraysize) - vfact(arraysize)*ufact(arraysize-1)) &
                          /(vfact(arraysize)*vfact(arraysize-1) + 1._dl)

       do i=1,(arraysize-1)
         j = arraysize - i
         yddot(j) = vfact(j)*yddot(j+1) + ufact(j)
       end do

!------splint.

       eind = 1
       do i=1,arraysize
         if (x.ge.xvalues(i)) eind = i
       end do

       if (eind.ge.arraysize) eind = arraysize - 1
       if (eind.lt.1) eind = 1

       parh = xvalues(eind+1) - xvalues(eind)
       para = (xvalues(eind+1) - x)/parh
       parb = (x - xvalues(eind))/parh

       ypred = para*yvalues(eind) + parb*yvalues(eind+1) &
               + ((para**3 - para)*yddot(eind) + (parb**3 - parb)*yddot(eind+1)) &
               *parh**2/6._dl

       deallocate(yddot)
       deallocate(ufact)
       deallocate(vfact)


       return


       end subroutine interpneff


!-----------------------------------------------------------

       subroutine cubic_spline(xvalues,yvalues,yddot,yp1,yp2)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Performs cubic spline using points.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: xvalues
       real(dl), dimension(:), intent(in) :: yvalues
       real(dl), dimension(:), intent(out) :: yddot
       real(dl), optional, intent(in) :: yp1
       real(dl), optional, intent(in) :: yp2


!------local variables.

       integer i,j !indicies

       integer arraysize !size of xvalues
       real(dl), dimension(:), allocatable :: ufact
       real(dl), dimension(:), allocatable :: vfact
       real(dl) wfact
       real(dl) sig1
       real(dl) sig2


!------procedure.


!------spline.

       arraysize = size(xvalues)

       allocate(ufact(arraysize)) !assign size of array ufact
       allocate(vfact(arraysize)) !assign size of array vfact

       yddot = 0._dl
       ufact = 0._dl
       vfact = 0._dl
       wfact = 0._dl

       if (present(yp1)) then
         vfact(1) = -0.5_dl
         ufact(1) = 3._dl/(xvalues(2) - xvalues(1)) &
                    *((yvalues(2) - yvalues(1))/(xvalues(2) - xvalues(1)) - yp1)
       else
         vfact(1) = 0._dl
         ufact(1) = 0._dl
       end if
       
       do i=2,(arraysize-1)
         sig1 = (xvalues(i) - xvalues(i-1))/(xvalues(i+1) - xvalues(i-1))
         sig2 = sig1*vfact(i-1) + 2._dl
         vfact(i) = (sig1 - 1._dl)/sig2
         wfact = &
              (yvalues(i+1) - yvalues(i))/(xvalues(i+1) - xvalues(i)) &
              - (yvalues(i) - yvalues(i-1))/(xvalues(i) - xvalues(i-1))
         ufact(i) = (6._dl*wfact/(xvalues(i+1) - xvalues(i-1)) &
              - sig1*ufact(i-1))/sig2
       end do

       if (present(yp2)) then
         vfact(arraysize) = 0.5_dl
         ufact(arraysize) = 3._dl/(xvalues(arraysize) - xvalues(arraysize-1)) &
                            *(yp2 - (yvalues(arraysize) - yvalues(arraysize-1)) &
                            /(xvalues(arraysize) - xvalues(arraysize-1)))
       else
         vfact(arraysize) = 0._dl
         ufact(arraysize) = 0._dl
       end if
       
       yddot(arraysize) = (ufact(arraysize) - vfact(arraysize)*ufact(arraysize-1)) &
                          /(vfact(arraysize)*vfact(arraysize-1) + 1._dl)

       do i=1,(arraysize-1)
         j = arraysize - i
         yddot(j) = vfact(j)*yddot(j+1) + ufact(j)
       end do

       deallocate(ufact)
       deallocate(vfact)


       return


       end subroutine cubic_spline


!-----------------------------------------------------------

       real(dl) function cubic_splint(xvalues,yvalues,yddot,x)


!------linkages.

!      called by - [subroutine] recom


!-----remarks.

!     Performs cubic splint using points.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: xvalues
       real(dl), dimension(:), intent(in) :: yvalues
       real(dl), dimension(:), intent(in) :: yddot
       real(dl), intent(in) :: x


!------local variables.

       integer i !index
       integer eind

       integer arraysize !size of xvalues
       real(dl) para
       real(dl) parb
       real(dl) parh


!------procedure.


!------splint.

       arraysize = size(xvalues)

       eind = 1
       do i=1,arraysize
         if (x.ge.xvalues(i)) eind = i
       end do

       if (eind.ge.arraysize) eind = arraysize - 1
       if (eind.lt.1) eind = 1

       parh = xvalues(eind+1) - xvalues(eind)
       para = (xvalues(eind+1) - x)/parh
       parb = (x - xvalues(eind))/parh

       cubic_splint = para*yvalues(eind) + parb*yvalues(eind+1) &
                      + ((para**3 - para)*yddot(eind) + (parb**3 - parb)*yddot(eind+1)) &
                      *parh**2/6._dl


       return


       end function cubic_splint


!-----------------------------------------------------------

       subroutine polint(xa,ya,x,y,dy)


!------linkages.

!      called by - [subroutine] scattering


!------remarks.

!      Performs polynomial interpolation using input points.
!      Maximum order: 9


       implicit none


!------throughput variables.

       real(kind = dl), dimension(:), intent(in) :: xa !input x array
       real(kind = dl), dimension(:), intent(in) :: ya !input y array
       real(kind = dl), intent(in) :: x !input x value
       real(kind = dl), intent(out) :: y !output y value
       real(kind = dl), intent(out) :: dy !error in y


!------local variables.

       integer i,m,ns,n !indicies
       real(kind = dl) dif, dift !differences between x and xa(i)
       real(kind = dl) ho, hp, w, den !quantities for c and d

       real(kind = dl), dimension(10) :: c !coefficients for interp
       real(kind = dl), dimension(10) :: d !coefficients for interp


!------procedure.

       n = size(xa)

       ns = 1
       dif = abs(x - xa(1))
       do i=1,n
         dift = abs(x - xa(i))
         if (dift.lt.dif) then
           ns = i
           dif = dift
         end if
         c(i) = ya(i)
         d(i) = ya(i)
       end do

       y = ya(ns)
       ns = ns - 1
       do m=1,n-1
         do i=1,n-m
           ho = xa(i) - x
           hp = xa(i+m) - x
           w = c(i+1) - d(i)
           den = ho - hp
           den = w/den
           c(i) = ho*den
           d(i) = hp*den
         end do
         if ((2*ns).lt.(n-m)) then
           dy = c(ns + 1)
         else
           dy = d(ns)
           ns = ns - 1
         end if
         y = y + dy
       end do


       return


       end subroutine polint


!-----------------------------------------------------------

       integer function indpolint(xa,x,ord)
       

!------linkages.

!      called by - [subroutine] scattering


!------remarks.

!      Finds the index needed for the interpolating arrays of polint.
!      Uses bisection method.


       implicit none


!------throughput variables.

       real(kind = dl), dimension(:), intent(in) :: xa !abscissas
       real(kind = dl), intent(in) :: x !input x value
       !integer, intent(out) :: j !output index for interpolating at x
       integer, intent(in) :: ord !order of interpolating polynomial


!------local variables.

       integer i, n !indices
       integer jl, jm, ju !indicies

       n = size(xa)

       jl = 0
       ju = n + 1

       do i=1,n
         if ((ju-jl).eq.1) exit
         jm = (ju + jl)/2
         if ((xa(n).ge.xa(1)).eqv.(x.ge.xa(jm))) then
           jl = jm
         else
           ju = jm
         end if
       end do !i

       if (jl.lt.(1 + ord/2)) then
         indpolint = 1
       else if (jl.gt.(n - 1 - ord/2)) then
         indpolint = n - ord
       else
         indpolint = jl - ord/2
       end if       


       return


       end function indpolint


       end module interp
