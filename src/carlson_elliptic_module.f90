!*******************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 2/7/2016
!
!  Carlson symmetric forms of elliptic integrals.
!
!  These routines are refactored versions of the ones from [SLATEC](http://www.netlib.org/slatec/).
!  They have been converted into modern Fortran, and the documentation has been
!  converted to FORD syntax.

    module carlson_elliptic_module

    use iso_fortran_env, only: error_unit, wp => real64

    implicit none

    private

    !**************************************************************
    !>
    !  Machine constants (replaces the old SLATEC [D1MACH](http://www.netlib.org/slatec/src/d1mach.f) function)
    !
    !  The traditional D1MACH constants are:
    !  * `D1MACH( 1) = B**(EMIN-1)`,           the smallest positive magnitude.
    !  * `D1MACH( 2) = B**EMAX*(1 - B**(-T))`, the largest magnitude.
    !  * `D1MACH( 3) = B**(-T)`,               the smallest relative spacing.
    !  * `D1MACH( 4) = B**(1-T)`,              the largest relative spacing.
    !  * `D1MACH( 5) = LOG10(B)`

    real(wp),dimension(5),parameter :: d1mach = &
        [  tiny(1.0_wp), &
           huge(1.0_wp), &
           real(radix(1.0_wp),wp)**(-digits(1.0_wp)), &
           epsilon(1.0_wp), &
           log10(real(radix(1.0_wp),wp)) ]

    !**************************************************************

    public :: drc,drd,drf,drj

    public :: elliptic_functions_test

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation of the Carlson elliptic integral:
!  $$ R_C(x,y) = \frac{1}{2} \int_{0}^{\infty} (t+x)^{-1/2}(t+y)^{-1} dt $$
!  where \(x\ge0\) and \(y>0\).
!
!  The duplication theorem is iterated until the variables are nearly equal,
!  and the function is then expanded in Taylor series to fifth order.
!  Logarithmic, inverse circular, and inverse hyperbolic functions can be
!  expressed in terms of DRC.
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### DRC special cases
!
!  $$
!    \begin{array}{rcll}
!     R_C(x,x+z) + R_C(y,y+z) &=& R_C(0,z)             & x>0, y>0, ~\mathrm{and}~ z>0 ~\mathrm{and}~ x y = z^2 \\
!     R_C(0,1/4)              &=& R_C(1/16,1/8) = \pi  & \\
!     R_C(9/4,2)              &=& \ln(2)               &
!    \end{array}
!  $$
!
!### Special functions via DRC
!
!  $$
!  \begin{array}{rll}
!     \ln(x)        &= (x-1) R_C \left( \left( \frac{1+x}{2} \right)^2, x \right) & x>0 \\
!     \sin^{-1}(x)  &=  x R_C ( (1-x)^2 ,1 ) & -1 \le x \le 1 \\
!     \cos^{-1}(x)  &= \sqrt{1-x^2} R_C(x^2,1)  & 0 \le x \le 1 \\
!     \tan^{-1}(x)  &= x R_C(1,1+x^2) & -\infty \lt x \lt \infty \\
!     \cot^{-1}(x)  &= R_C(x^2, x^2+1) & 0 \le x \lt \infty \\
!     \sinh^{-1}(x) &= x R_C(1+x^2, 1) & -\infty \lt x \lt \infty \\
!     \cosh^{-1}(x) &= \sqrt{x^2-1} R_C(x^2,1)  & x \ge 1 \\
!     \tanh^{-1}(x) &= x R_C(1,1-x^2) & -1 < x < 1 \\
!     \coth^{-1}(x) &= R_C(x^2,x^2-1) & x > 1 \\
!  \end{array}
!  $$
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drc.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(wp) function drc(x,y,ier)

    implicit none

    real(wp),intent(in) :: x    !! nonnegative variable
    real(wp),intent(in) :: y    !! positive variable
    integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! `IER = 0`: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! `IER > 0`: Abnormal termination of the routine:
                                !! `IER = 1`: `x<0 or y<=0`
                                !! `IER = 2`: `x+y<LOLIM`
                                !! `IER = 3`: `max(x,y) > UPLIM`

    character(len=16) :: xern3 , xern4 , xern5
    real(wp) :: lamda, mu , s , sn , xn , yn

    real(wp),parameter :: errtol = (d1mach(3)/16.0_wp)**(1.0_wp/6.0_wp)
        !! Determines the accuracy of the answer.
        !!
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! `16 * ERRTOL ** 6 / (1 - 2 * ERRTOL)`.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 2.0e-17),
        !! (3.0e-3, 2.0e-14),
        !! (1.0e-2, 2.0e-11),
        !! (3.0e-2, 2.0e-8),
        !! (1.0e-1, 2.0e-5)
        !!
        !! The accuracy of the computed approximation to the inte-
        !! gral can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the trunca-
        !! tion error there will be round-off error, but in prac-
        !! tice the total error from both sources is usually less
        !! than the amount given in the table.
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.
    real(wp),parameter :: lolim  = 5.0_wp*d1mach(1) !! Lower limit of valid arguments
    real(wp),parameter :: uplim  = d1mach(2)/5.0_wp !! Upper limit of valid arguments
    real(wp),parameter :: c1     = 1.0_wp/7.0_wp
    real(wp),parameter :: c2     = 9.0_wp/22.0_wp

    !initialize:
    drc = 0.0_wp

    ! check for errors:
    if ( x<0.0_wp .or. y<=0.0_wp ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write(error_unit,'(a)') &
            'drc: x<0 .or. y<=0 where x = '//xern3//' and y = '//xern4
        return
    endif

    if ( max(x,y)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') uplim
        write(error_unit,'(a)') &
            'drc: max(x,y)>uplim where x = '//&
            xern3//' y = '//xern4//' and uplim = '//xern5
        return
    endif

    if ( x+y<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') lolim
        write(error_unit,'(a)') &
            'drc: x+y<lolim where x = '//xern3//&
            ' y = '//xern4//' and lolim = '//xern5
        return
    endif

    ier = 0
    xn = x
    yn = y

    do
        mu = (xn+yn+yn)/3.0_wp
        sn = (yn+mu)/mu - 2.0_wp
        if ( abs(sn)<errtol ) exit
        lamda = 2.0_wp*sqrt(xn)*sqrt(yn) + yn
        xn = (xn+lamda)*0.250_wp
        yn = (yn+lamda)*0.250_wp
    end do

    s = sn*sn*(0.30_wp+sn*(c1+sn*(0.3750_wp+sn*c2)))
    drc = (1.0_wp+s)/sqrt(mu)

    end function drc
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 2nd kind:
!  $$ R_D(x,y,z) = \frac{3}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-3/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(x+y>0\), and \(z>0\).
!
!  If \(x=0\) or \(y=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRD Special Comments
!
!  $$
!    \begin{array}{rl}
!      R_D(x,y,z) + R_D(y,z,x) + R_D(z,x,y) = \frac{3}{\sqrt{x y z}}  &  x>0, y>0, z>0
!    \end{array}
!  $$
!
!### Special functions via DRD and DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      E(\phi,k) = \sin \phi  R_F(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!      -\frac{k^2}{3} \sin^3 \phi R_D(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!    $$
!    When \( \phi = \pi /2 \) the integral is complete:
!    $$
!      \begin{array}{rcl}
!       E(k) &=& R_F(0,1-k^2 ,1) - \frac{k^2}{3} R_D(0,1-k^2 ,1) \\
!            &=& \int_{0}^{\pi/2} \sqrt{1-k^2 \sin^2 \phi}  d \phi
!      \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \mathrm{EL2}(x,k_c,a,b) = ax R_F(1,1+k_c^2 x^2 ,1+x^2 )
!        + \frac{1}{3}(b-a) x^3 R_D(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Legendre form of alternative ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \begin{array}{rcl}
!      D(q,k) &=& \int_{0}^{q} \sin^2 p (1-k^2 \sin^2 p)^{-1/2} dp \\
!             &=& \frac{1}{3} (\sin^3 q) R_D(\cos^2 q,1-k^2 \sin^2 q,1)
!      \end{array}
!    $$
!
!  * Lemniscate constant B:
!
!    $$
!      \begin{array}{rcl}
!      B &=& \int_{0}^{1} s^2 (1-s^4)^{-1/2} ds \\
!        &=& \frac{1}{3} R_D (0,2,1)
!      \end{array}
!    $$
!
!  * Heuman's LAMBDA function:
!
!    $$
!      \begin{array}{rcl}
!      \frac{\pi}{2} \Lambda_0(a,b) &=&
!      \sin b \left(R_F(0,\cos^2 a,1)-\frac{1}{3} \sin^2 a
!      R_D(0,\cos^2 a,1) \right) R_F(\cos^2 b,1-\cos^2 a \sin^2 b,1) \\
!      & &-\frac{1}{3} \cos^2 a \sin^3 b R_F(0,\cos^2 a,1)
!      R_D(\cos^2 b,1-\cos^2 a \sin^2 b,1)
!      \end{array}
!    $$
!
!  * Jacobi ZETA function:
!
!    $$
!      \begin{array}{rcl}
!      Z(b,k) &=& \frac{k^2}{3} \sin b R_F(\cos^2 b,1-k^2 \sin^2 b,1)
!               R_D(0,1-k^2 ,1)/R_F(0,1-k^2 ,1) \\
!                & & -\frac{k^2}{3} \sin^3 b R_D(\cos^2 b,1-k^2 \sin^2 b,1)
!      \end{array}
!    $$
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Modify calls to XERMSG to put in standard form.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drd.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(wp) function drd(x,y,z,ier)

    implicit none

    real(wp),intent(in) :: x    !! nonnegative variable (\(x+y>0\))
    real(wp),intent(in) :: y    !! nonnegative variable (\(x+y>0\))
    real(wp),intent(in) :: z    !! positive variable
    integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! `IER = 0`: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! `IER > 0`: Abnormal termination of the routine:
                                !! `IER = 1`: `min(x,y) < 0`
                                !! `IER = 2`: `min(x + y, z ) < LOLIM`
                                !! `IER = 3`: `max(x,y,z) > UPLIM`

    character(len=16) :: xern3 , xern4 , xern5 , xern6
    real(wp) :: epslon, ea , eb , ec , ed , ef , lamda
    real(wp) :: mu , power4 , sigma , s1 , s2 , xn , xndev
    real(wp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

    real(wp),parameter :: errtol = (d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
        !! Determines the accuracy of the answer.
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! `3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2`.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 4.0e-18),
        !! (3.0e-3, 3.0e-15),
        !! (1.0e-2, 4.0e-12),
        !! (3.0e-2, 3.0e-9),
        !! (1.0e-1, 4.0e-6)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    real(wp),parameter :: lolim  = 2.0_wp/(d1mach(2))**(2.0_wp/3.0_wp) !! Lower limit of valid arguments
    real(wp),parameter :: tuplim = (0.10_wp*errtol)**(1.0_wp/3.0_wp)/&
                                    d1mach(1)**(1.0_wp/3.0_wp)
    real(wp),parameter :: uplim  = tuplim**2  !! Upper limit of valid arguments
    real(wp),parameter :: c1     = 3.0_wp/14.0_wp
    real(wp),parameter :: c2     = 1.0_wp/6.0_wp
    real(wp),parameter :: c3     = 9.0_wp/22.0_wp
    real(wp),parameter :: c4     = 3.0_wp/26.0_wp

    ! initialize:
    drd = 0.0_wp

    ! check for errors:
    if ( min(x,y)<0.0_wp ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write(error_unit,'(a)') 'drd: min(x,y)<0 where x = '//xern3// &
                                ' and y = '//xern4
        return
    endif

    if ( max(x,y,z)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drd: max(x,y,z)>uplim where x = '// &
                                xern3//' y = '//xern4//' z = '//xern5// &
                                ' and uplim = '//xern6
        return
    endif

    if ( min(x+y,z)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drd: min(x+y,z)<lolim where x = '// &
                                xern3//' y = '//xern4//' z = '//xern5// &
                                ' and lolim = '//xern6
        return
    endif

    ier    = 0
    xn     = x
    yn     = y
    zn     = z
    sigma  = 0.0_wp
    power4 = 1.0_wp

    do
        mu     = (xn+yn+3.0_wp*zn)*0.20_wp
        xndev  = (mu-xn)/mu
        yndev  = (mu-yn)/mu
        zndev  = (mu-zn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        sigma  = sigma + power4/(znroot*(zn+lamda))
        power4 = power4*0.250_wp
        xn     = (xn+lamda)*0.250_wp
        yn     = (yn+lamda)*0.250_wp
        zn     = (zn+lamda)*0.250_wp
    end do

    ea  = xndev*yndev
    eb  = zndev*zndev
    ec  = ea - eb
    ed  = ea - 6.0_wp*eb
    ef  = ed + ec + ec
    s1  = ed*(-c1+0.250_wp*c3*ed-1.50_wp*c4*zndev*ef)
    s2  = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
    drd = 3.0_wp*sigma + power4*(1.0_wp+s1+s2)/(mu*sqrt(mu))

    end function drd
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 1st kind:
!  $$ R_F(x,y,z) = \frac{1}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-1/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(z\ge0\), and at most one of
!  them is \(=0\).
!
!  If \(x=0\), \(y=0\), or \(z=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRF Special Comments
!
!  $$
!    \begin{array}{rl}
!    R_F(x,x+z,x+w) + R_F(y,y+z,y+w) = R_F(0,z,w)
!    & x>0, y>0, z>0, x y = z w
!    \end{array}
!  $$
!
!### Special functions via DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!    \begin{array}{rl}
!      F(\phi,k) &= \sin \phi R_F( \cos^2 \phi,1-k^2 \sin^2 \phi,1) \\
!           K(k) &= R_F(0,1-k^2 ,1) = \int_{0}^{\pi/2} (1-k^2 sin^2 \phi )^{-1/2} d \phi
!    \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!      \mathrm{EL1}(x,k_c) = x R_F(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Lemniscate constant A:
!
!    $$
!      A = \int_{0}^{1} (1-s^4 )^{-1/2}    ds = R_F(0,1,2) = R_F(0,2,1)
!    $$
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC))
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drf.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(wp) function drf(x,y,z,ier)

    implicit none

    real(wp),intent(in) :: x    !! nonnegative variable
    real(wp),intent(in) :: y    !! nonnegative variable
    real(wp),intent(in) :: z    !! nonnegative variable
    integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! `IER = 0`: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! `IER > 0`: Abnormal termination of the routine:
                                !! `IER = 1`: `min(x,y,z) < 0`
                                !! `IER = 2`:` min(x+y,x+z,y+z) < LOLIM`
                                !! `IER = 3`: `max(x,y,z) > UPLIM`

    character(len=16) :: xern3 , xern4 , xern5 , xern6
    real(wp) :: epslon, e2 , e3 , lamda
    real(wp) :: mu , s , xn , xndev
    real(wp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

    real(wp),parameter :: errtol = (4.0_wp*d1mach(3))**(1.0_wp/6.0_wp)
        !! Determines the accuracy of the answer.
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! `ERRTOL ** 6 / (4 * (1-ERRTOL)`.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 3.0e-19),
        !! (3.0e-3, 2.0e-16),
        !! (1.0e-2, 3.0e-13),
        !! (3.0e-2, 2.0e-10),
        !! (1.0e-1, 3.0e-7)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    real(wp),parameter :: lolim  = 5.0_wp*d1mach(1) !! Lower limit of valid arguments
    real(wp),parameter :: uplim  = d1mach(2)/5.0_wp !! Upper limit of valid arguments
    real(wp),parameter :: c1     = 1.0_wp/24.0_wp
    real(wp),parameter :: c2     = 3.0_wp/44.0_wp
    real(wp),parameter :: c3     = 1.0_wp/14.0_wp

    ! initialize:
    drf = 0.0_wp

    ! check for errors:
    if ( min(x,y,z)<0.0_wp ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write(error_unit,'(a)') 'drf: min(x,y,z)<0 where x = '// &
                xern3//' y = '//xern4//' and z = '//xern5
        return
    endif

    if ( max(x,y,z)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drf: max(x,y,z)>uplim where x = '// &
                xern3//' y = '//xern4//' z = '//xern5// &
                ' and uplim = '//xern6
        return
    endif

    if ( min(x+y,x+z,y+z)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drf: min(x+y,x+z,y+z)<lolim where x = '//xern3// &
                ' y = '//xern4//' z = '//xern5//' and lolim = '//xern6
        return
    endif

    ier = 0
    xn  = x
    yn  = y
    zn  = z

    do
        mu     = (xn+yn+zn)/3.0_wp
        xndev  = 2.0_wp - (mu+xn)/mu
        yndev  = 2.0_wp - (mu+yn)/mu
        zndev  = 2.0_wp - (mu+zn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        xn     = (xn+lamda)*0.250_wp
        yn     = (yn+lamda)*0.250_wp
        zn     = (zn+lamda)*0.250_wp
    end do
    e2  = xndev*yndev - zndev*zndev
    e3  = xndev*yndev*zndev
    s   = 1.0_wp + (c1*e2-0.10_wp-c2*e3)*e2 + c3*e3
    drf = s/sqrt(mu)

    end function drf
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 3rd kind:
!  $$
!    R_J(x,y,z,p) = \frac{3}{2} \int_{0}^{\infty}
!                   (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} (t+p)^{-1} dt
!  $$
!  where \(x\ge0\), \(y\ge0\), and \(z\ge0\),
!  and at most one of them \(=0\), and \(p>0\).
!
!  If \(x=0\) or \(y=0\) or \(z=0\), then the
!  integral is COMPLETE.
!
!  The duplication theorem is iterated
!  until the variables are nearly equal, and the function is
!  then expanded in Taylor series to fifth order.
!
!### DRJ Special Comments
!
!  $$
!  R_J(x,x+z,x+w,x+p) + R_J(y,y+z,y+w,y+p) +
!    (a-b) R_J(a,b,b,a) + \frac{3}{\sqrt{a}} = R_J(0,z,w,p)
!  $$
!  where:
!  $$
!      x>0, y>0, z>0, w>0, p>0
!  $$
!  and:
!  $$
!     \begin{array}{rl}
!         x y &= z w           \\
!           a &= p^2 (x+y+z+w) \\
!           b &= p (p+x) (p+y) \\
!       b - a &= p (p-z) (p-w)
!     \end{array}
!  $$
!  The sum of the third and
!  fourth terms on the left side is \( 3 R_C(a,b) \) .
!
!### Special functions via DRJ and DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 3rd kind:
!
!  $$
!      \begin{array}{rcl}
!      P(\phi,k,n) &=& \int_{0}^{\phi} (1+n \sin^2 \theta )^{-1}
!                    (1-k^2 \sin^2 \theta )^{-1/2} d \theta \\
!                  &=& \sin \phi R_F(\cos^2 \phi, 1-k^2 sin^2 \phi,1)
!                    -\frac{n}{3} \sin^3 \phi R_J(\cos^2 \phi,1-k^2 \sin^2 \phi, 1,1+n \sin^2 \phi)
!      \end{array}
!  $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind:
!
!  $$
!   \begin{array}{rcl}
!      \mathrm{EL3}(x,k_c,p) &=& x R_F(1,1+k_c^2 x^2 ,1+x^2 ) +
!      \frac{1}{3}(1-p) x^3 R_J(1,1+k_c^2 x^2 ,1+x^2 ,1+p x^2 ) \\
!      \mathrm{CEL}(k_c,p,a,b) &=& a R_F(0,k_c^2 ,1) +
!      \frac{1}{3}(b-pa) R_J(0,k_c^2 ,1,p)
!   \end{array}
!  $$
!
!  * Heuman's LAMBDA function:
!
!  $$
!   \begin{array}{rcl}
!   L(a,b,p) &=&  ( \cos^2 a \sin b \cos b /(1-\cos^2 a \sin^2 b )^{1/2} )
!          (\sin p R_F(\cos^2 p ,1-\sin^2 a \sin^2 p ,1) \\
!          & ~ & +
!          (\sin^2 a \sin^3 p /(3(1-\cos^2 a \sin^2 b )))
!          R_J(\cos^2 p ,1-\sin^2 a \sin^2 p ,1,1-
!          \sin^2 a \sin^2 p /(1-\cos^2 a \sin^2 b ))) \\
!   \frac{\pi}{2} \Lambda_0(a,b) &=& L(a,b,\pi/2) \\
!         &=& \cos^2 a  \sin b \cos b (1-\cos^2 a \sin^2 b )^{-1/2}
!         R_F(0,\cos^2 a ,1) \\
!         &~& + (1/3) \sin^2 a \cos^2 a
!         \sin b \cos b (1-\cos^2 a \sin^2 b )^{-3/2}
!         R_J(0,\cos^2 a ,1,\cos^2 a \cos^2 b /(1-\cos^2 a \sin^2 b ))
!   \end{array}
!  $$
!
!  * Jacobi ZETA function:
!
!  $$
!      Z(b,k) = \frac{k^2}{3} \sin b \cos b (1-k^2 \sin^2 b)^{1/2}
!               \frac{R_J(0,1-k^2 ,1,1-k^2 \sin^2 b)}{R_F (0,1-k^2 ,1)}
!  $$
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC)).
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drj.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(wp) function drj(x,y,z,p,ier)

    implicit none

    real(wp),intent(in) :: x    !! nonnegative variable
    real(wp),intent(in) :: y    !! nonnegative variable
    real(wp),intent(in) :: z    !! nonnegative variable
    real(wp),intent(in) :: p    !! positive variable
    integer,intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! `IER = 0`: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! `IER = 1`: `min(x,y,z) < 0.0_wp`
                                !! `IER = 2`: `min(x+y,x+z,y+z,p) < LOLIM`
                                !! `IER = 3`: `max(x,y,z,p) > UPLIM`

    character(len=16) xern3 , xern4 , xern5 , xern6 , xern7
    real(wp) :: alfa , beta , ea , eb , ec , e2 , e3, epslon
    real(wp) :: lamda , mu , pn , pndev
    real(wp) :: power4 , sigma , s1 , s2 , s3 , xn , xndev
    real(wp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

    real(wp),parameter :: errtol = (d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
        !! Determines the accuracy of the answer
        !!
        !! the value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! "machine precision".
        !!
        !! Relative error due to truncation of the series for DRJ
        !! is less than `3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2`.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than}
        !! (1.0e-3, 4.0e-18),
        !! (3.0e-3, 3.0e-15),
        !! (1.0e-2, 4.0e-12),
        !! (3.0e-2, 3.0e-9),
        !! (1.0e-1, 4.0e-6)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    ! LOLIM and UPLIM determine the valid range of X, Y, Z, and P
    real(wp),parameter :: lolim  = (5.0_wp*d1mach(1))**(1.0_wp/3.0_wp)
        !! not less than the cube root of the value
        !! of LOLIM used in the routine for [[DRC]].
    real(wp),parameter :: uplim  = 0.30_wp*(d1mach(2)/5.0_wp)**(1.0_wp/3.0_wp)
        !! not greater than 0.3 times the cube root of
        !! the value of UPLIM used in the routine for DRC.

    real(wp),parameter :: c1     = 3.0_wp/14.0_wp
    real(wp),parameter :: c2     = 1.0_wp/3.0_wp
    real(wp),parameter :: c3     = 3.0_wp/22.0_wp
    real(wp),parameter :: c4     = 3.0_wp/26.0_wp

    !initialize:
    drj = 0.0_wp

    ! check for errors:
    if ( min(x,y,z)<0.0_wp ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write(error_unit,'(a)') 'drj: min(x,y,z)<0 where x = '// &
            xern3//' y = '//xern4//' and z = '//xern5
        return
    endif

    if ( max(x,y,z,p)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') p
        write (xern7,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drj: max(x,y,z,p)>uplim where x = '// &
            xern3//' y = '//xern4//' z = '//xern5//' p = '// &
            xern6//' and uplim = '//xern7
        return
    endif

    if ( min(x+y,x+z,y+z,p)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') p
        write (xern7,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drj: min(x+y,x+z,y+z,p)<lolim where x = '//xern3// &
            ' y = '//xern4//' z = '//xern5//' p = '//xern6// &
            ' and lolim = '//xern7
        return
    endif

    ier    = 0
    xn     = x
    yn     = y
    zn     = z
    pn     = p
    sigma  = 0.0_wp
    power4 = 1.0_wp

    do
        mu     = (xn+yn+zn+pn+pn)*0.20_wp
        xndev  = (mu-xn)/mu
        yndev  = (mu-yn)/mu
        zndev  = (mu-zn)/mu
        pndev  = (mu-pn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev),abs(pndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        alfa   = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
        alfa   = alfa*alfa
        beta   = pn*(pn+lamda)*(pn+lamda)
        sigma  = sigma + power4*drc(alfa,beta,ier)
        power4 = power4*0.250_wp
        xn     = (xn+lamda)*0.250_wp
        yn     = (yn+lamda)*0.250_wp
        zn     = (zn+lamda)*0.250_wp
        pn     = (pn+lamda)*0.250_wp
    end do

    ea  = xndev*(yndev+zndev) + yndev*zndev
    eb  = xndev*yndev*zndev
    ec  = pndev*pndev
    e2  = ea - 3.0_wp*ec
    e3  = eb + 2.0_wp*pndev*(ea-ec)
    s1  = 1.0_wp + e2*(-c1+0.750_wp*c3*e2-1.50_wp*c4*e3)
    s2  = eb*(0.50_wp*c2+pndev*(-c3-c3+pndev*c4))
    s3  = pndev*ea*(c2-pndev*c3) - c2*pndev*ec
    drj = 3.0_wp*sigma + power4*(s1+s2+s3)/(mu*sqrt(mu))

    end function drj
!*******************************************************************************

!*******************************************************************************
!>
!  Unit tests for the [[carlson_elliptic_module]] module.

    subroutine elliptic_functions_test()

    implicit none

    real(wp) :: x,y,z,p,w
    real(wp) :: d,truth
    integer :: ier
    real(wp),parameter :: pi = acos(-1.0_wp)

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' elliptic_functions_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) 'DRC ERRTOL=',(d1mach(3)/16.0_wp)**(1.0_wp/6.0_wp)
    write(*,*) 'DRD ERRTOL=',(d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
    write(*,*) 'DRF ERRTOL=',(4.0_wp*d1mach(3))**(1.0_wp/6.0_wp)
    write(*,*) 'DRJ ERRTOL=',(d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
    write(*,*) ''

    !DRC:

    x = 8.0_wp
    y = 2.0_wp
    z = 4.0_wp
    d = drc(x,x+z,ier) + drc(y,y+z,ier)
    truth = drc(0.0_wp,z,ier)
    write(*,'(A60,E20.9)') 'DRC(X,X+Z) + DRC(Y,Y+Z) = DRC(0,Z) error = ', d-truth

    x = 0.0_wp
    y = 1.0_wp / 4.0_wp
    d = drc(x,y,ier)
    truth = pi
    write(*,'(A60,E20.9)') 'DRC(0,1/4) error = ', d-truth

    x = 1.0_wp / 16.0_wp
    y = 1.0_wp / 8.0_wp
    d = drc(x,y,ier)
    truth = pi
    write(*,'(A60,E20.9)') 'DRC(1/16,1/8) error = ', d-truth

    x = 9.0_wp / 4.0_wp
    y = 2.0_wp
    d = drc(x,y,ier)
    truth = log(2.0_wp)
    write(*,'(A60,E20.9)') 'DRC(9/4,2) error = ', d-truth

    !DRD:

    x = 10.0_wp
    y = 20.0_wp
    z = 30.0_wp
    d = drd(x,y,z,ier) + drd(y,z,x,ier) + drd(z,x,y,ier)
    truth = 3.0_wp / sqrt(x * y * z)
    write(*,'(A60,E20.9)') 'DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y) error = ', d-truth

    !DRF:

    x = 10.0_wp
    y = 20.0_wp
    z = 100.0_wp
    w = 2.0_wp
    d = DRF(X,X+Z,X+W,ier) + DRF(Y,Y+Z,Y+W,ier)
    truth = DRF(0.0_wp,Z,W,ier)
    write(*,'(A60,E20.9)') 'DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W) error = ', d-truth

    !DRJ:

    !...

    write(*,*) ''

    end subroutine elliptic_functions_test
!*******************************************************************************

!*******************************************************************************
    end module carlson_elliptic_module
!*******************************************************************************
