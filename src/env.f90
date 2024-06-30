module env
    use iso_fortran_env, only: wp => real64
    implicit none

    real(wp), parameter :: epsilon = tiny(1.0_wp)
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
end module env