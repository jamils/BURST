
       module transvar


!------linkages.

!      used by - [program] main


!-----remarks.

!     Contains a list of variables and related functions for main


       use mainvar


       implicit none


       real(dl) nuratiotol !tolerance for \nu ratio quantities.

       real(dl), dimension(:,:,:), allocatable :: woccfrac
       real(dl), dimension(:,:,:), allocatable :: prec_trans

       real(dl) :: gfhbar = 0.20668190937_dl !G_F**2/hbar (in s^{-1}MeV^{-5})
       real(dl) :: sin2thw = 0.23_dl

       integer :: ordpolint = 5


       type(bin_scheme) glgndr !for Gauss-Legendre Integration
       type(bin_scheme) glagur !for Gauss-Laguerre Integration


       !generic arrays used in calculating rates:
       real(dl), dimension(:,:,:,:,:), allocatable :: ratescatt
       real(dl), dimension(:,:,:,:,:), allocatable :: frsratescatt
       real(dl), dimension(:,:,:,:), allocatable :: rateannih
       real(dl), dimension(:,:,:,:), allocatable :: frsrateannih

       real(dl), dimension(:,:,:,:,:), allocatable :: scatt
       real(dl), dimension(:,:,:,:,:), allocatable :: frsscatt
       real(dl), dimension(:,:,:,:), allocatable :: annih
       real(dl), dimension(:,:,:,:), allocatable :: frsannih

       real(dl), dimension(:), allocatable :: glegai
       real(dl), dimension(:), allocatable :: glegwi

       real(dl), dimension(:), allocatable :: glegai2
       real(dl), dimension(:), allocatable :: glegwi2

       real(dl), dimension(:), allocatable :: glagai2
       real(dl), dimension(:), allocatable :: glagwi2

       real(dl), dimension(:), allocatable :: glegai3
       real(dl), dimension(:), allocatable :: glegwi3

       real(dl), dimension(:), allocatable :: e2a
       real(dl), dimension(:), allocatable :: e2w

       real(dl), dimension(:), allocatable :: e3a
       real(dl), dimension(:), allocatable :: e3w

       real(dl), dimension(:), allocatable :: glagai3
       real(dl), dimension(:), allocatable :: glagwi3

       real(dl), dimension(:), allocatable :: ggenai2
       real(dl), dimension(:), allocatable :: ggenwi2

       real(dl), dimension(:), allocatable :: e3r2a
       real(dl), dimension(:), allocatable :: e3r2w

       real(dl), dimension(:), allocatable :: algndrout
       real(dl), dimension(:), allocatable :: wlgndrout

       real(dl), dimension(:), allocatable :: alagurout
       real(dl), dimension(:), allocatable :: wlagurout

       real(dl), dimension(:), allocatable :: algndrin
       real(dl), dimension(:), allocatable :: wlgndrin

       real(dl), dimension(:), allocatable :: alagurin
       real(dl), dimension(:), allocatable :: wlagurin

       real(dl), dimension(:), allocatable :: agenout
       real(dl), dimension(:), allocatable :: wgenout

       real(dl), dimension(:), allocatable :: agenin
       real(dl), dimension(:), allocatable :: wgenin

       real(dl), dimension(:), allocatable :: eouta
       real(dl), dimension(:), allocatable :: eoutw

#ifdef prllel

       real(dl), dimension(:,:,:,:,:), allocatable :: rsprllel
       real(dl), dimension(:,:,:,:,:), allocatable :: frsrsprllel
       real(dl), dimension(:,:,:,:), allocatable :: raprllel
       real(dl), dimension(:,:,:,:), allocatable :: frsraprllel
       real(dl), dimension(:,:,:), allocatable :: dfdtprllel

#endif


       contains


!-----------------------------------------------------------

       subroutine interp_occfrac(x,y,xval,yval)


!------linkages.

!      called by - [program] scattering_nunu
!      calls     - none

!------remarks.

!      Subroutine to interpolate occupation probabilities.


       use interp


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: x !x values
       real(dl), dimension(:,:,:), intent(in) :: y !y values
       real(dl), intent(in) :: xval !input x value
       real(dl), dimension(:,:), intent(out) :: yval !computed y values


!------local variables.

       integer n,m !indicies
       integer j
       real(kind = dl) errory !error used in evaluating polint
       real(dl) xbar, amp, sigmax


!------procedure.


!------interpolated values.

       if (xval.lt.0._dl) then
#ifdef prllel
         write (*,*) xval, rank
#else
         write (*,*) xval
#endif
         stop
       end if

       if (xval.le.2._dl*x(size(x))) then
       !if (xval.le.40._dl) then
         j = indpolint(x,xval,ordpolint)
         !j = int(xval/x(size(x))*real(nbins, kind=dl)) - 1
         !j = max(j,1)
         !j = min(j,size(x)-ordpolint)
         do n=1,size(yval, 2)
           do m=1,size(yval, 1)
             call polint(x(j:j+ordpolint),y(m,n,j:j+ordpolint),xval,yval(m,n),errory)
             !06Oct2015 EG: insert safety statement for MSW alterations
             if ((isnan(yval(m,n))).or.(yval(m,n).lt.-300._dl) &
                 .or.(yval(m,n).gt.0._dl)) then
               yval(m,n) = -300._dl
             end if
             !write (*,*)
             !x(j:j+ordpolint),y(m,n,j:j+ordpolint),xval,yval(m,n),errory
           end do !m
         end do !n
       else if ((xval.gt.2._dl*x(size(x))).and.(xval.lt.300._dl)) then
         do n=1,size(yval, 2)
           do m=1,size(yval, 1)
             call polint(x(size(x)-1:),y(m,n,size(x)-1:),xval,yval(m,n),errory)
             !write (*,*)
             !x(size(x)-1:),y(m,n,size(x)-1:),xval,yval(m,n),errory
           end do !m
         end do !n
       else
         yval(:,:) = -300._dl
       end if


       !xbar = x(16)
       !amp = 1.0_dl/(exp(xbar) + 1._dl)
       !sigmax = x(2) - x(1)
       !yval(:,:) = log(1._dl/(exp(xval) + 1._dl) + amp*exp(-(xval-xbar)**2/2._dl/sigmax**2))


       return


       end subroutine interp_occfrac


       end module transvar
