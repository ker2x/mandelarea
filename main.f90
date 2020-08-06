program mandelarea
!http://idris.fr/formations/fortran
use aplot
implicit none

type(aplot_t)::p
real, dimension(:), allocatable :: x !x-axis
real, dimension(:), allocatable :: y    !y-axis

real, dimension(:), allocatable :: inside    !notinmset counter
real, dimension(:), allocatable :: outside   !inmset counter
!real, dimension(10)::y = (/ 2.0, 4.0, 6.0, 8.0, 6.0, 4.0, 2.0, 0.0, 2.0, 4.0 /)

real :: r, i    !real & imaginary to create a random complex
complex :: z    !the random complex to be initialized
integer :: b, j ! loop stuff

integer, parameter :: startiter = 10             !start at iteration startiter
integer, parameter :: itermultiplicator = 100     !multiply by itermultiplicator at each batch
integer, parameter :: batch = 2000                !number of batch to compute
integer, parameter :: loopmax = 200000            !montecarlo loop
integer :: currentiter = 0

print *, "Starting MandelArea"

p = initialize_plot()

!le code par ici
ALLOCATE(x(batch)) 
ALLOCATE(y(batch)) 
ALLOCATE(inside(batch)) 
ALLOCATE(outside(batch)) 

x = 0
y = 0
inside = 0
outside = 0 
currentiter =  startiter

DO b = 1, batch
    x(b) = currentiter
    print *, b, currentiter
    do j = 1, loopmax
        call random_number(r)
        call random_number(i)    

        z = CMPLX(r*2.5-2, i*2.6-1.3)    
        if(isInMSet(z, currentiter)) then
            inside(b) = inside(b) + 1
        else
            outside(b) = outside(b) + 1
        end if
    end do
    currentiter = currentiter + itermultiplicator
end do


    call add_dataset(p, x, outside / inside)
    !call set_yscale (p, minval(outside / inside), maxval(outside / inside))
    call set_xscale(p,real(startiter), real(currentiter))
    call set_yscale(p,2.5,4.0)
    !call set_xlogarithmic (p,2.0)
    call set_xlabel(p, "iteration")
    call set_ylabel(p, "ratio out/in")
    call set_title(p, "Mandelbrot area")
    call set_serieslabel(p, 0, "Mandelbrot area")
    call set_seriestype(p, 0, APLOT_STYLE_DOT)
    call display_plot(p)
    call destroy_plot(p)

DEALLOCATE(x)
DEALLOCATE(y)
DEALLOCATE(inside)
DEALLOCATE(outside)

CONTAINS 

pure function isInMSet(c, maxiter)
    complex, intent(in) :: c
    integer, intent(in) :: maxiter
    integer :: n
    complex :: z
    logical :: isInMSet
    
    z = CMPLX(0,0)
    n = 0
    IF(((ABS(c - CMPLX(-1,0))) < 0.25) .OR. ((ABS(1.0 - SQRT(1-(4*c)))) < 1.0)) THEN
        isInMSet = .TRUE.
    END IF

    DO WHILE(ABS(z) < 4 .AND. (n < maxiter))
        z = z*z+c
        n = n+1
    END do

    if(n >= maxiter) then
        isInMSet = .TRUE.
    else
        isInMset = .FALSE.
    end if    
end function isInMSet

end program mandelarea
