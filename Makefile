#
# Automagically generated by Approximatrix Simply Fortran 2.41
#
FC="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\gfortran.exe"
CC="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\gcc.exe"
AR="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\ar.exe"
WRC="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\windres.exe"
RM=rm -f

IDIR=

LDIR=


OPTFLAGS= -O3 -fgraphite-identity -floop-interchange -floop-strip-mine -floop-block -floop-parallelize-all -mtune=bdver4

SPECIALFLAGS=$(IDIR)

RCFLAGS=-O coff

PRJ_FFLAGS=-ffast-math -floop-parallelize-all -ftree-parallelize-loops=4 -fcheck=all -fopenmp -std=f2008

PRJ_CFLAGS=

PRJ_LFLAGS=-laplot -lappgraphics -lgdi32 -lcomdlg32 -lcomctl32 -luuid -loleaut32 -lole32 -lstdc++

FFLAGS=$(SPECIALFLAGS) $(OPTFLAGS) $(PRJ_FFLAGS) -Jmodules 

CFLAGS=$(SPECIALFLAGS) $(OPTFLAGS) $(PRJ_CFLAGS)

"build\main.o": ".\main.f90"
	@echo Compiling .\main.f90
	@$(FC) -c -o "build\main.o" $(FFLAGS) ".\main.f90"

clean: .SYMBOLIC
	@echo Deleting build\main.o and related files
	@$(RM) "build\main.o"
	@echo Deleting program.exe
	@$(RM) "program.exe"

"program.exe":  "build\main.o" "build\mandelarea.prj.target"
	@echo Generating program.exe
	@$(FC) -o "program.exe" -static -fopenmp "build\main.o" $(LDIR) $(PRJ_LFLAGS)

all: "program.exe" .SYMBOLIC

