program mandelarea

use aplot

implicit none

type(aplot_t)::p

real, dimension(:), allocatable :: x !x-axis
real, dimension(:), allocatable :: inside    !notinmset counter
real, dimension(:), allocatable :: outside   !inmset counter

real    :: r, i !real & imaginary to create a random complex
complex :: z    !the random complex to be initialized
integer :: b, j ! loop stuff

integer :: currentiter = 0
integer, parameter :: startiter = 100            !start at iteration startiter
integer, parameter :: iterstep = 100           !add(or mult by) iterstep at each batch
integer, parameter :: batch = 100               !number of point on x-axis
integer, parameter :: loopmax = 500000          !montecarlo loop

ALLOCATE(x(batch)) 
ALLOCATE(inside(batch)) 
ALLOCATE(outside(batch)) 

x = 0
inside = 0
outside = 0 
currentiter =  startiter

DO b = 1, batch
    x(b) = real(currentiter)
    !try some omp stuff in the inner loop ?
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(inside, outside,b,currentiter) private(z,r,i,x)
    DO j = 1, loopmax
        call random_number(r)
        call random_number(i)    
        !z = CMPLX(r*2.5-2, i*2.6-1.3)    
        z = CMPLX(r*4.0-2.0, i*4.0-2.0)    

        do while(abs(z) .GT. 2.0)
            call random_number(r)
            call random_number(i)    
            z = CMPLX(r*4.0-2.0, i*4.0-2.0)    
        end do

        if(isInMSet(z, currentiter)) then
            !$omp atomic update
            inside(b) = inside(b) + 1
        else
            !$omp atomic update
            outside(b) = outside(b) + 1
        end if
    end do
    !$OMP END PARALLEL DO
    print *, b, currentiter, inside(b), outside(b), (inside(b) / outside(b)) * 100.0
    currentiter = currentiter + iterstep
end do

p = initialize_plot()
call add_dataset(p, x, (inside / outside) * 100.0)
!call set_yscale(p,9.0,12.0)
!call set_xscale(p,real(startiter), real(currentiter))
!call set_xlogarithmic (p,2.0)
call set_xlabel(p, "iteration")
call set_ylabel(p, "area ratio")
call set_title(p, "Mandelbrot area")
call set_serieslabel(p, 0, "Mandelbrot area")
call set_seriestype(p, 0, APLOT_STYLE_DOT)
call display_plot(p)
call destroy_plot(p)

DEALLOCATE(x)
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
