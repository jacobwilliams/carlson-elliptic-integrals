!*******************************************************************************
!>
!  Unit tests for the [[carlson_elliptic_module]] module.

program carlson_elliptic_test

    use carlson_elliptic_module, wp => carlson_elliptic_module_wp

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

    ! write(*,*) 'DRC ERRTOL=',(d1mach(3)/16.0_wp)**(1.0_wp/6.0_wp)
    ! write(*,*) 'DRD ERRTOL=',(d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
    ! write(*,*) 'DRF ERRTOL=',(4.0_wp*d1mach(3))**(1.0_wp/6.0_wp)
    ! write(*,*) 'DRJ ERRTOL=',(d1mach(3)/3.0_wp)**(1.0_wp/6.0_wp)
    ! write(*,*) ''

    !DRC:

    x = 8.0_wp
    y = 2.0_wp
    z = 4.0_wp
    d = drc(x,x+z,ier) + drc(y,y+z,ier)
    truth = drc(0.0_wp,z,ier)
    write(*,'(A60,E20.9)') 'DRC(X,X+Z) + DRC(Y,Y+Z) = DRC(0,Z) error = ', d-truth
    if (ier /= 0) error stop ier

    x = 0.0_wp
    y = 1.0_wp / 4.0_wp
    d = drc(x,y,ier)
    truth = pi
    write(*,'(A60,E20.9)') 'DRC(0,1/4) error = ', d-truth
    if (ier /= 0) error stop ier

    x = 1.0_wp / 16.0_wp
    y = 1.0_wp / 8.0_wp
    d = drc(x,y,ier)
    truth = pi
    write(*,'(A60,E20.9)') 'DRC(1/16,1/8) error = ', d-truth
    if (ier /= 0) error stop ier

    x = 9.0_wp / 4.0_wp
    y = 2.0_wp
    d = drc(x,y,ier)
    truth = log(2.0_wp)
    write(*,'(A60,E20.9)') 'DRC(9/4,2) error = ', d-truth
    if (ier /= 0) error stop ier

    !DRD:

    x = 10.0_wp
    y = 20.0_wp
    z = 30.0_wp
    d = drd(x,y,z,ier) + drd(y,z,x,ier) + drd(z,x,y,ier)
    truth = 3.0_wp / sqrt(x * y * z)
    write(*,'(A60,E20.9)') 'DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y) error = ', d-truth
    if (ier /= 0) error stop ier

    !DRF:

    x = 10.0_wp
    y = 20.0_wp
    z = 100.0_wp
    w = 2.0_wp
    d = DRF(X,X+Z,X+W,ier) + DRF(Y,Y+Z,Y+W,ier)
    truth = DRF(0.0_wp,Z,W,ier)
    write(*,'(A60,E20.9)') 'DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W) error = ', d-truth
    if (ier /= 0) error stop ier

    !DRJ:
    d = drj(0.5_wp,0.5_wp,0.5_wp,2.0_wp,ier)
    truth = 1.11836068453037130354466230219139_wp
    write(*,'(A60,E20.9)') 'drj(0.5,0.5,0.5,2) error = ', d-truth
    if (ier /= 0) error stop ier

end program carlson_elliptic_test
