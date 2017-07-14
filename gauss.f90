
       module gauss


!------linkages.

!      used by - [program] main


!------remarks.

!      Contains relevant subroutines and functions from gglr.f90
!      used to calculate Gauss quadrature abscissas and weights.


       use mainvar


       implicit none


       contains


!-------------------------------------------------------------------     

       subroutine cgqf (gqr,t,wts,gaussflag,alpha,beta)


!------linkages.

!      called by - [program] main
!      calls     - [subroutine] cdgqf


!------remarks.

!      cgqf computes knots and weights of a Gauss quadrature formula.


!      Discussion:

!      Only simple knots are produced.
!      The following rules are available:

!      1, Legendre,             (-1,1)       1.0
!      2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!      3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!      4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!      5, Generalized Laguerre, (0,inf)     (x-a)^alpha*exp(-b*(x-a))
!      6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!      7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!      8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      16 February 2010


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:
!
!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       integer, intent(in) :: gqr !Gauss quadrature rule.
       real(dl), dimension(:), intent(out) :: t !knots.
       real(dl), dimension(:), intent(out) :: wts !weights.
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.


!------procedure.


       gaussflag = .true.

!      Compute the Gauss quadrature formula for default endpoint values of a and b.

       if (present(alpha)) then
         if (present(beta)) then
           call cdgqf (gqr,t,wts,gaussflag,alpha,beta)
         else
           call cdgqf (gqr,t,wts,gaussflag,alpha)
         end if
       else
         call cdgqf (gqr,t,wts,gaussflag)
       end if


       return


       end subroutine cgqf


!-------------------------------------------------------------------     

       subroutine cdgqf(gqr,t,wts,gaussflag,alpha,beta)


!------linkages.

!      called by - [subroutine] cgqf
!      calls     - [subroutine] parchk, class_matrix, sgqf


!------remarks.

!      cdgqf computes a Gauss quadrature formula with default (a,b) and simple knots.


!      Discussion:

!      This routine computes all the knots and weights of a Gauss quadrature
!      formula with a classical weight function with default values for a and b,
!      and only simple knots.  The following rules are available:

!      1, Legendre,             (a,b)       1.0
!      2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!      3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!      4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!      5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!      6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!      7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!      8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      04 January 2010


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.

!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       integer, intent(in) :: gqr !Gauss quadrature rule.
       real(dl), dimension(:), intent(out) :: t !knots.
       real(dl), dimension(:), intent(out) :: wts !weights.
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.

       integer nt !size of array t
       real(dl), dimension(:), allocatable :: aj
       real(dl), dimension(:), allocatable :: bj
       real(dl) zemu


!------procedure.


       nt = size(t)

       allocate(aj(nt)) !assign dimension of aj
       allocate(bj(nt)) !assign dimension of bj


!      Check input parameters for specific rule.
       if (present(alpha)) then
         if (present(beta)) then
           call parchk (gqr,2*nt,gaussflag,alpha,beta)
         else
           call parchk (gqr,2*nt,gaussflag,alpha)
         end if
       else
         call parchk (gqr,2*nt,gaussflag)
       end if
       if (.not.gaussflag) return


!      Get the Jacobi matrix and zero-th moment.
       if (present(alpha)) then
         if (present(beta)) then
           call class_matrix (gqr,aj,bj,zemu,gaussflag,alpha,beta)
         else
           call class_matrix (gqr,aj,bj,zemu,gaussflag,alpha)
         end if
       else
         call class_matrix (gqr,aj,bj,zemu,gaussflag)
       end if
       if (.not.gaussflag) return


!      Compute the knots and weights.
       call sgqf(aj,bj,zemu,t,wts,gaussflag)


       deallocate(aj)
       deallocate(bj)


       return


       end subroutine cdgqf


!-------------------------------------------------------------------     

       subroutine parchk (gqr,m,gaussflag,alpha,beta)


!------linkages.

!      called by - [subroutine] cdgqf,scqf
!      calls     - [subroutine] parchk


!------remarks.

!      parchk checks parameters ALPHA and BETA for classical weight functions. 


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      27 December 2009


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       integer, intent(in) :: gqr !rule.
       integer, intent(in) :: m !order of highest moment (only needed if gqr = 8).
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.

       real(dl) tmp


!------procedure.


       if (gqr.le.0) gaussflag = .false.

!      Check alpha for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.

       if ((present(alpha)).and.(gqr.ge.3)) then
         if (alpha.le.-1._dl) gaussflag = .false.
       end if

!      Check beta for Jacobi.

       if ((present(beta)).and.(gqr.eq.4)) then
         if (beta.le.-1._dl) gaussflag = .false.
       end if

!      Check alpha and beta for rational.

       if ((present(alpha)).and.(present(beta)).and.(gqr.eq.8)) then
         tmp = alpha + beta + m + 1._dl
         if ((tmp.ge.0._dl).or.(tmp.le.beta)) then
           gaussflag = .false.
         end if
       end if


       return


       end subroutine parchk


!-------------------------------------------------------------------     

       subroutine class_matrix (gqr,aj,bj,zemu,gaussflag,alpha,beta)


!------linkages.

!      called by - [subroutine] cdgqf,scqf
!      calls     - [subroutine] parchk
!                  [function]   r8_gamma


!------remarks.

!      class_matrix computes the Jacobi matrix for a quadrature rule.


!      Discussion:

!      This routine computes the diagonal aj and sub-diagonal bj
!      elements of the order M tridiagonal symmetric Jacobi matrix
!      associated with the polynomials orthogonal with respect to
!      the weight function specified by gqr.

!      For weight functions 1-7, m elements are defined in bj even
!      though only m-1 are needed.  For weight function 8, bj(m) is
!      set to zero.

!      The zero-th moment of the weight function is returned in zemu.

!      The following rules are available:

!      1, Legendre,             (a,b)       1.0
!      2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!      3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!      4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!      5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!      6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!      7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!      8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      27 December 2009


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.



       implicit none


!------throughput variables.

       integer, intent(in) :: gqr !Gauss quadrature rule.
       real(dl), dimension(:), intent(out) :: aj !the diagonal of the Jacobi matrix.
       real(dl), dimension(:), intent(out) :: bj !the subdiagonal of the Jacobi matrix.
       real(dl), intent(out) :: zemu !zero-th moment.
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.

       integer m !order of the Jacobi matrix.
       real(dl) a2b2
       real(dl) ab
       real(dl) aba
       real(dl) abi
       real(dl) abj
       real(dl) abti
       real(dl) apone
       integer i !index
       real(dl) ireal
       real(dl) temp
       real(dl) temp2


!------procedure.


       !temp = epsilon(temp)
       temp = epsilon(1._dl)

       m = size(aj)

       if (present(alpha)) then
         if (present(beta)) then
           call parchk(gqr,2*m-1,gaussflag,alpha,beta)
         else
           call parchk(gqr,2*m-1,gaussflag,alpha)
         end if
       else
           call parchk(gqr,2*m-1,gaussflag)
       end if
       if (.not.gaussflag) return


       temp2 = 0.5_dl

       if (500._dl*temp.lt.abs((r8_gamma(temp2))**2 - pi)) then
         gaussflag = .false.
         return
       end if


       if (gqr.eq.1) then

         ab = 0._dl
         zemu = 2._dl/(ab + 1._dl)
         aj = 0._dl
         do i=1,m
           ireal = real(i, kind=dl)
           abi = ireal + ab*real(mod(i,2), kind=dl)
           abj = 2._dl*ireal + ab
           bj(i) = abi**2/(abj**2 - 1._dl)
         end do
         bj(:) = sqrt(bj(:))

       else if (gqr.eq.2) then

         zemu = pi
         aj = 0._dl
         bj(:) = 0.5_dl
         bj(1) =  sqrt(0.5_dl)

       else if (gqr.eq.3) then

         ab = 2._dl*alpha
         zemu = 2.0_dl**(ab + 1._dl)*r8_gamma(alpha + 1._dl)**2 &
                /r8_gamma(ab + 2._dl)
         aj = 0._dl
         bj(1) = 1._dl/(2._dl*alpha + 3._dl)
         do i=2,m
           ireal = real(i, kind=dl)
           bj(i) = ireal*(ireal + ab)/(4._dl*(ireal + alpha)**2 - 1._dl)
         end do
         bj(:) = sqrt(bj(:))

       else if (gqr.eq.4) then

         ab = alpha + beta
         abi = 2._dl + ab
         zemu = 2._dl**(ab + 1._dl)*r8_gamma(alpha + 1._dl) &
                *r8_gamma(beta + 1._dl)/r8_gamma(abi)
         aj(1) = (beta - alpha)/abi
         bj(1) = 4._dl*(1._dl + alpha)*(1._dl + beta) &
                 /((abi + 1._dl)*abi**2)
         a2b2 = beta**2 - alpha**2
         do i=2,m
           ireal = real(i, kind=dl)
           abi = 2._dl*ireal + ab
           aj(i) = a2b2/((abi - 2._dl)*abi)
           abi = abi**2
           bj(i) = 4._dl*ireal*(ireal + alpha)*(ireal + beta)*(ireal + ab) &
                   /((abi - 1._dl)*abi)
         end do
         bj(:) = sqrt(bj(:))

       else if (gqr.eq.5) then

         zemu = r8_gamma(alpha + 1._dl)
         do i=1,m
           ireal = real(i, kind=dl)
           aj(i) = 2._dl*ireal - 1._dl + alpha
           bj(i) = ireal*(ireal + alpha)
         end do
         bj(:) =  sqrt(bj(:))

       else if (gqr.eq.6) then

         zemu = r8_gamma((alpha + 1._dl)/2._dl)
         aj = 0._dl
         do i=1,m
           bj(i) = (real(i, kind=dl) + alpha*real(mod(i,2), kind=dl))/2._dl
         end do
         bj(:) = sqrt(bj(:))

       else if (gqr.eq.7) then

         ab = alpha
         zemu = 2._dl/(ab + 1._dl)
         aj = 0._dl
         do i=1,m
           ireal = real(i, kind=dl)
           abi = ireal + ab*real(mod(i,2), kind=dl)
           abj = 2*ireal + ab
           bj(i) = abi**2/(abj**2 - 1._dl)
         end do
         bj(:) = sqrt(bj(:))

       else if (gqr.eq.8) then

         ab = alpha + beta
         zemu = r8_gamma(alpha + 1._dl)*r8_gamma(-(ab + 1._dl)) &
               /r8_gamma(-beta)
         apone = alpha + 1._dl
         aba = ab*apone
         aj(1) = -apone/(ab + 2._dl)
         bj(1) = -aj(1)*(beta + 1._dl)/(ab + 2._dl)/(ab + 3._dl)
         do i=2,m
           ireal = real(i, kind=dl)
           abti = ab + 2._dl*ireal
           aj(i) = aba + 2._dl*(ab + ireal)*(ireal - 1._dl)
           aj(i) = -aj(i)/abti/(abti - 2._dl)
         end do
         do i=2,(m-1)
           ireal = real(i, kind=dl)
           abti = ab + 2._dl*ireal
           bj(i) = ireal*(alpha + ireal)/(abti - 1._dl)*(beta + ireal) &
                   /abti**2*(ab + ireal)/(abti + 1._dl)
         end do
         bj(m) = 0._dl
         bj(:) = sqrt(bj(:))

       end if


       return


       end subroutine class_matrix


!-------------------------------------------------------------------     


       subroutine sgqf(aj,bj,zemu,t,wts,gaussflag)


!------linkages.

!      called by - [subroutine] cdgqf
!      calls     - [subroutine] imtqlx


!------remarks.

!      sgqf computes knots and weights of a Gauss Quadrature formula.


!      Discussion:

!      This routine computes all the knots and weights of a Gauss quadrature
!      formula with simple knots from the Jacobi matrix and the zero-th
!      moment of the weight function, using the Golub-Welsch technique.


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      04 January 2010


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: aj !diagonal of the Jacobi matrix.
       real(dl), dimension(:), intent(inout) :: bj !subdiagonal of the Jacobi matrix
       real(dl), intent(in) :: zemu !zeroth-order moment
       real(dl), dimension(:), intent(inout) :: t !knots.
       real(dl), dimension(:), intent(inout) :: wts !weights.
       logical, intent(inout) :: gaussflag !Flag for failure.


!------local variables.

       integer i !index


!------procedure.


!      Exit if the zero-th moment is not positive.

       if (zemu.le.0._dl) then
         gaussflag = .false.
         return
       end if


!      Set up vectors for IMTQLX.

       t(:) = aj(:)
       wts = 0._dl
       wts(1) = sqrt(zemu)


!      Diagonalize the Jacobi matrix.

       call imtqlx(t,bj,wts,gaussflag)

       wts(:) = wts(:)**2


       return


       end subroutine sgqf


!-------------------------------------------------------------------     


       subroutine imtqlx(d,e,z,gaussflag)


!------linkages.

!      called by - [subroutine] sgqf
!      calls     - none


!------remarks.

!      imtqlx diagonalizes a symmetric tridiagonal matrix.

!      Discussion:

!      This routine is a slightly modified version of the EISPACK routine to 
!      perform the implicit QL algorithm on a symmetric tridiagonal matrix. 

!      The authors thank the authors of EISPACK for permission to use this
!      routine. 

!      It has been modified to produce the product Q' * Z, where Z is an input 
!      vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!      The changes consist (essentially) of applying the orthogonal 
!      transformations directly to Z as they are generated.


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      27 December 2009


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.

!      Roger Martin, James Wilkinson,
!      The Implicit QL Algorithm,
!      Numerische Mathematik,
!      Volume 12, Number 5, December 1968, pages 377-383.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(inout) :: d !diagonal entries of the matrix.
       real(dl), dimension(:), intent(out) :: e !subdiagonal entries of the matrix.
       real(dl), dimension(:), intent(inout) :: z !vector.
       logical, intent(inout) :: gaussflag !Flag for failure.


!------local variables.

        integer n !size of array d
        real(dl) b
        real(dl) c
        real(dl) f
        real(dl) g
        integer i
        integer ii
        integer, parameter :: itn = 30
        integer j
        integer k
        integer l
        integer m
        integer mml
        real(dl) p
        real(dl) prec
        real(dl) r
        real(dl) s


!------procedure.


       n = size(d)

       prec = epsilon(prec)

       if (n.eq.1) return

       e(n) = 0._dl


       do l=1,n
         j = 0
         do
           do m=l,n
             if (m.eq.n) exit
             if (abs(e(m)).le.prec*(abs(d(m)) + abs(d(m+1)))) exit
           end do
           p = d(l)
           if (m.eq.l) exit
           if (itn.le.j) then
             gaussflag = .false.
             return
           end if
           j = j + 1
           g = (d(l+1) - p)/(2._dl*e(l))
           r = sqrt(g**2 + 1._dl)
           g = d(m) - p + e(l)/(g + sign(r,g))
           s = 1._dl
           c = 1._dl
           p = 0._dl
           mml = m - l
           do ii=1,mml
             i = m - ii
             f = s*e(i)
             b = c*e(i)
             if (abs(g).le.abs(f)) then
               c = g/f
               r = sqrt(c**2 + 1._dl)
               e(i+1) = f*r
               s = 1._dl/r
               c = c*s
             else
               s = f/g
               r = sqrt(s**2 + 1._dl)
               e(i+1) = g*r
               c = 1._dl/r
               s = s*c
             end if
             g = d(i+1) - p
             r = (d(i) - g)*s + 2._dl*c*b
             p = s*r
             d(i+1) = g+p
             g = c*r - b
             f = z(i+1)
             z(i+1) = s*z(i) + c*f
             z(i) = c*z(i) - s*f
           end do
           d(l) = d(l) - p
           e(l) = g
           e(m) = 0._dl
         end do
       end do


!      Sorting:
       do ii=2,n
         i = ii - 1
         k = i
         p = d(i)
         do j=ii,n
           if (d(j).lt.p) then
             k = j
             p = d(j)
           end if
         end do
         if (k.ne.i) then
           d(k) = d(i)
           d(i) = p
           p = z(i)
           z(i) = z(k)
           z(k) = p
         end if
       end do


       return


       end subroutine imtqlx


!-------------------------------------------------------------------     

       subroutine cscgqf(gqr,t,wts,a,b,gaussflag,alpha,beta)


!------linkages.

!      called by - [program] main
!      calls     - [subroutine] scqf


!------remarks.

!      cscgqf computes knots and weights of a scaled Gauss quadrature formula.

!      Discussion:

!      The user may specify the interval (a,b).  Only simple knots are produced.
!      The following rules are available:

!      1, Legendre,             (a,b)       1.0
!      2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!      3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!      4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!      5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!      6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!      7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!      8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      16 February 2010


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.

!      Reference:
!
!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       integer, intent(in) :: gqr !Gauss quadrature rule.
       real(dl), dimension(:), intent(out) :: t !knots.
       real(dl), dimension(:), intent(out) :: wts !weights.
       real(dl), intent(in) :: a !interval starting point
       real(dl), intent(in) :: b !interval ending point
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.

       integer i !index
       integer nt !size of array t
       integer, dimension(:), allocatable :: mlt
       integer, dimension(:), allocatable :: ndx


!------procedure.


       gaussflag = .true.

       nt = size(t)

!      Prepare to scale the quadrature formula to other weight function with 
!      valid a and b:

       allocate(mlt(nt)) !assign dimension of mlt
       mlt = 1

       allocate(ndx(nt)) !assign dimension of ndx
       do i = 1, nt 
         ndx(i) = i
       end do

!      scale the quadrature formula to other weight function:
       if (present(alpha)) then
         if (present(beta)) then
           call scqf(t,mlt,wts,ndx,wts,t,gqr,a,b,gaussflag,alpha,beta)
         else
           call scqf(t,mlt,wts,ndx,wts,t,gqr,a,b,gaussflag,alpha)
         end if
       else
         call scqf(t,mlt,wts,ndx,wts,t,gqr,a,b,gaussflag)
       end if


       deallocate(mlt)
       deallocate(ndx)


       return


       end subroutine cscgqf


!-------------------------------------------------------------------     

       subroutine scqf(t,mlt,wts,ndx,swts,st,gqr,a,b,gaussflag,alpha,beta)


!------linkages.

!      called by - [subroutine] cdgqf
!      calls     - [subroutine] parchk


!------remarks.

!      scqf scales a quadrature formula to a nonstandard interval.

!      Discussion:

!      The arrays WTS and SWTS may coincide.  The arrays T and ST may coincide.
!      The following rules are available:

!      1, Legendre,             (a,b)       1.0
!      2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!      3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!      4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!      5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!      6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!      7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!      8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta


!      Licensing:

!      This code is distributed under the GNU LGPL license. 


!      Modified:

!      27 December 2009


!      Author:

!      Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      Sylvan Elhay, Jaroslav Kautsky,
!      Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!      Interpolatory Quadrature,
!      ACM Transactions on Mathematical Software,
!      Volume 13, Number 4, December 1987, pages 399-415.


       implicit none


!------throughput variables.

       real(dl), dimension(:), intent(in) :: t !original knots.
       integer, dimension(:), intent(in) :: mlt !multiplicity of the knots.
       real(dl), dimension(:), intent(in) :: wts !WTS(NWTS), the weights.
       integer, dimension(:), intent(in) :: ndx !used to index the array wts.  
       real(dl), dimension(:), intent(out) :: swts !SWTS(NWTS), the scaled weights.
       real(dl), dimension(:), intent(out) :: st !scaled knots.
       integer, intent(in) :: gqr !rule.
       real(dl), intent(in) :: a !interval starting point
       real(dl), intent(in) :: b !interval ending point
       logical, intent(inout) :: gaussflag !Flag for failure.
       real(dl), optional, intent(in) :: alpha !the value of alpha, if needed.
       real(dl), optional, intent(in) :: beta !the value of beta, if needed.


!------local variables.

       real(dl) al
       real(dl) be
       integer i,k,l !indicies
       real(dl) p
       real(dl) shft
       real(dl) slp
       real(dl) temp
       real(dl) tmp


!------procedure.


       !temp = epsilon(temp)
       temp = epsilon(1._dl)


       if (present(alpha)) then
         if (present(beta)) then
           call parchk (gqr,1,gaussflag,alpha,beta)
         else
           call parchk (gqr,1,gaussflag,alpha)
         end if
       else
         call parchk (gqr,1,gaussflag)
       end if
       if (.not.gaussflag) return

 
       if (gqr.eq.1) then

         if (abs(b-a).le.temp) then
           gaussflag = .false.
         else
           al = 0._dl
           be = 0._dl
           shft = (a + b)/2._dl
           slp = (b - a)/2._dl
         end if

       else if (gqr.eq.2) then

         if (abs(b - a).le.temp) then
           gaussflag = .false.
         else
           al = -0.5_dl
           be = -0.5_dl
           shft = (a + b)/2._dl
           slp = (b - a)/2._dl
         end if

       else if (gqr.eq.3) then

         if (abs(b - a).le.temp) then
           gaussflag = .false.
         else
           al = alpha
           be = alpha
           shft = (a + b)/2._dl
           slp = (b - a)/2._dl
         end if

       else if (gqr.eq.4) then

         if (abs(b - a).le.temp) then
           gaussflag = .false.
         else
           al = alpha
           be = beta
           shft = (a + b)/2._dl
           slp = (b - a)/2._dl
         end if

       else if (gqr.eq.5) then

         if (b.le.0._dl) then
           gaussflag = .false.
         else
           al = alpha
           be = 0._dl
           shft = a
           slp = 1._dl/b
         end if

       else if (gqr.eq.6) then

         if (b.le.0._dl) then
           gaussflag = .false.
         else
           al = alpha
           be = 0._dl
           shft = a
           slp = 1._dl/sqrt(b)
         end if

       else if (gqr.eq.7) then

         if (abs(b - a).le.temp) then
           gaussflag = .false.
         else
           al = alpha
           be = 0._dl
           shft = (a + b)/2._dl
           slp = (b - a)/2._dl
         end if

       else if (gqr.eq.8) then

         if ((a + b).le.0._dl) then
           gaussflag = .false.
         else
           al = alpha
           be = beta
           shft = a
           slp = a + b
         end if

       end if
       if (.not.gaussflag) return


       p = slp**(al + be + 1._dl)
       do k =1,size(t)
         st(k) = shft + slp*t(k)
         l = abs(ndx(k))
         if (l.ne.0) then
           tmp = p
           do i=l,(l+mlt(k)-1)
             swts(i) = wts(i)*tmp
             tmp = tmp*slp
           end do
         end if
       end do


       return


       end subroutine scqf


!-------------------------------------------------------------------     


       real(dl) function r8_gamma(x)


!------linkages.

!      called by - [subroutine] class_matrix
!      calls     - none


!------remarks.

!      r8_gamma evaluates Gamma(x) for a real argument.

!      Discussion:

!      This routine calculates the gamma function for a real argument x.

!      Computation is based on an algorithm outlined in reference 1.
!      The program uses rational functions that approximate the gamma
!      function to at least 20 significant decimal digits.  Coefficients
!      for the approximation over the interval (1,2) are unpublished.
!      Those for the approximation for 12 <= x are from reference 2.


!      Modified:

!      11 February 2008


!      Author:

!      Original FORTRAN77 version by William Cody, Laura Stoltz.
!      FORTRAN90 version by John Burkardt.


!      Reference:

!      William Cody,
!      An Overview of Software Development for Special Functions,
!      in Numerical Analysis Dundee, 1975,
!      edited by GA Watson,
!      Lecture Notes in Mathematics 506,
!      Springer, 1976.

!      John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!      Charles Mesztenyi, John Rice, Henry Thatcher,
!      Christoph Witzgall,
!      Computer Approximations,
!      Wiley, 1968,
!      LC: QA297.C64.


       implicit none


!------throughput variables.

       real(dl), intent(in) :: x !argument of the function.


!------local variables.

       real(dl), dimension(7) :: c = (/ &
                                      -1.910444077728e-03_dl, &
                                      8.4171387781295e-04_dl, &
                                      -5.952379913043012e-04_dl, &
                                      7.93650793500350248e-04_dl, &
                                      -2.777777777777681622553e-03_dl, &
                                      8.333333333333333331554247e-02_dl, &
                                      5.7083835261e-03_dl/)
       real(dl), dimension (8) :: p = (/ &
                                       -1.71618513886549492533811_dl, &
                                       2.47656508055759199108314e+01_dl, &
                                       -3.79804256470945635097577e+02_dl, &
                                       6.29331155312818442661052e+02_dl, &
                                       8.66966202790413211295064e+02_dl, &
                                       -3.14512729688483675254357e+04_dl, &
                                       -3.61444134186911729807069e+04_dl, &
                                       6.64561438202405440627855e+04_dl/)
       real(dl), dimension (8) :: q = (/ &
                                       -3.08402300119738975254353e+01_dl, &
                                       3.15350626979604161529144e+02_dl, &
                                       -1.01515636749021914166146e+03_dl, &
                                       -3.10777167157231109440444e+03_dl, &
                                       2.25381184209801510330112e+04_dl, &
                                       4.75584627752788110767815e+03_dl, &
                                       -1.34659959864969306392456e+05_dl, &
                                       -1.15132259675553483497211e+05_dl/)
       real(dl), parameter :: eps = 2.22e-16_dl
       real(dl) fact
       integer i
       integer n
       logical parity
       real(dl) res
       real(dl), parameter :: sqrtpi = 0.9189385332046727417803297_dl
       real(dl) sum
       real(dl), parameter :: xbig = 171.624_dl
       real(dl) xden
       real(dl), parameter :: xinf = 1.e+30_dl
       real(dl), parameter :: xminin = 2.23e-308_dl
       real(dl) xnum
       real(dl) y
       real(dl) y1
       real(dl) ysq
       real(dl) z


!------procedure.


       parity = .false.
       fact = 1._dl
       n = 0
       y = x

!      Argument is negative:
       if (y.le.0._dl) then
         y = -x
         y1 = aint(y)
         res = y - y1
         if (res.ne.0._dl) then
           if (y1.ne.aint(y1*0.5_dl)*2._dl) then
             parity = .true.
           end if
           fact = -pi/sin(pi*res)
           y = y + 1._dl
         else
           res = xinf
           r8_gamma = res
           return
         end if
       end if


!      Argument is positive:
       if (y.lt.eps) then
!        Argument < EPS:
         if (xminin.le.y) then
           res = 1._dl/y
         else
           res = xinf
           r8_gamma = res
           return
         end if
       else if (y.lt.12._dl) then
         y1 = y
!        0.0 < argument < 1.0:
         if (y.lt.1._dl) then
           z = y
           y = y + 1._dl
!          1.0 < argument < 12.0:
!          Reduce argument if necessary:
         else
           n = int(y) - 1
           y = y - real(n, kind=dl)
           z = y - 1._dl
         end if
!        Evaluate approximation for 1.0 < argument < 2.0:
         xnum = 0._dl
         xden = 1._dl
         do i=1,8
           xnum = (xnum + p(i))*z
           xden = xden*z + q(i)
         end do
         res = xnum/xden + 1._dl
!        Adjust result for case  0.0 < argument < 1.0:
         if (y1.lt.y) then
           res = res/y1
!        Adjust result for case 2.0 < argument < 12.0:
         else if (y.lt.y1) then
           do i=1,n
             res = res*y
             y = y + 1._dl
           end do
         end if
       else
!      Evaluate for 12.0 <= argument:
         if (y.le.xbig) then
           ysq = y**2
           sum = c(7)
           do i=1,6
             sum = sum/ysq + c(i)
           end do
           sum = sum/y - y + sqrtpi
           sum = sum + (y - 0.5_dl)*log(y)
           res = exp(sum)
         else
           res = xinf
           r8_gamma = res
           return
         end if
       end if


!      Final adjustments and return:
       if (parity) res = -res
       if (fact.ne.1._dl) res = fact/res
       r8_gamma = res


       return


       end function r8_gamma


       end module gauss
