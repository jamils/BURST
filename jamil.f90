module jamil
    
    implicit none

    contains

    REAL(KIND=8) function cvc(x)

        REAL(KIND=8) :: x

        cvc = -0.000128145 + (8.66854 * (10.**(-6)) * x) + &
        (8.37184 * (10.**(-6)) * x**2) + (5.32329 * (10.**(-8)) * x**3) + &
        (1.62786 * (10.**(-10)) * x**4) - (2.43836 * (10.**(-13)) * x**5) + &
        (1.43109 * (10.**(-16)) * x**6)

    end

end module