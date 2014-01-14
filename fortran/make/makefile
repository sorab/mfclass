#must use gfortran 4.8.2 or later to compile this program.

# Makefile example for compiling MODFLOW-2005:

# Makefile for non-USGS distribution
#
PROGRAM = mf2015
SOURCEDIR=../src
VPATH=${SOURCEDIR}
.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

# Define the Fortran compile flags

F90FLAGS=
F90= gfortran

# 
# Define the C compile flags
# -D_UF defines UNIX naming conventions for mixed language compilation.
# 
CFLAGS= -D_UF -O3
CC= gcc


# Define the libraries

SYSLIBS= -lc
USRLIB  = 

# Define all Fortran object files which make up MODFLOW

OBJECTS = \
        utl7u1.o \
        parutl7.o \
        loose_modules.o \
        xmdlib_2.o \
        sms_modules.o \
        xmd_modules.o \
        package.o \
        pcgu_modules.o \
        sparse.o \
        package_wel.o \
        package_ghb.o \
        model.o \
        cross.o \
        soln.o \
        sms7.o \
        xmd.o \
        pcgu7.o \
        mf2015.o \
        


# Define Task Function

all: mf2015

# Define what mf2015

mf2015: $(OBJECTS)
	-$(F90) $(F90FLAGS) -o mf2015 $(OBJECTS) $(USRLIB) $(SYSLIBS)

# Object codes

.f.o:
	$(F90) $(F90FLAGS) -c $<
	
.f90.o:
	$(F90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<
	
#  end
