#======================================================================================================================
# Change these as needed

# Are we making a static library, a shared library, or an executable?
#TO_BUILD := STATICLIB
TO_BUILD := SHAREDLIB
#TO_BUILD := EXE

# Is this a debug or a release build of CShore?
BUILD_VERSION := DEBUG
#BUILD_VERSION := RELEASE
#BUILD_VERSION := VALGRIND

# Specify the CShore input and output approach:
# input and output via arguments to/from calling program (ARGINOUT)
# input and output via arguments to/from calling program, also file output for checking purposes (ARGINBOTHOUT)
# input and output via ASCII files read by CShore executable or by calling program (FILEINOUT)
#INPUT_AND_OUTPUT := ARGINOUT
#INPUT_AND_OUTPUT := ARGINBOTHOUT
INPUT_AND_OUTPUT := FILEINOUT
#======================================================================================================================

# Fortran directives
FC := gfortran
AR := ar
ARFLAGS := cr
LD := gfortran
LDFLAGS :=

MAKEFILE := Makefile
HEADER := cshore.h
STATICLIB := libcshore.a
SHAREDLIB := libcshore.so
EXE := cshore

SRC := $(wildcard *.f03)
OBJS := $(SRC:.f03=.o)
MODSRC := cshore_USACE.f03
MOD := CshoreShared.mod

# DFM Got some not-very-useful (?) "Conversion from REAL(4) to REAL(8)" warnings with -Wconversion-extra, so not using

ifeq ($(BUILD_VERSION), DEBUG)
	FFLAGS := -std=f2008 -g -DDEBUG -O0 -Wall -Wextra -Wno-compare-reals -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std -Wsurprising -Waliasing -Wshadow -Wuninitialized -Wcast-align -Wnull-dereference -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal -fcheck=all -ffree-line-length-none -fno-underscoring -cpp -I.
	CPEXT := .DEBUG.LINUX

else ifeq ($(BUILD_VERSION), RELEASE)
	FFLAGS := -std=f2008 -O3 -ffree-line-length-none -fno-underscoring -cpp -I.
	CPEXT := .RELEASE.LINUX

else ifeq ($(BUILD_VERSION), VALGRIND)
	FFLAGS := -pg -std=f2008 -g -O0 -Wall -Wextra -Wno-compare-reals -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std -Wsurprising -Waliasing -Wshadow -Wuninitialized -Wcast-align -Wnull-dereference -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal -fcheck=all -ffree-line-length-none -fno-underscoring -cpp -I.
	LDFLAGS := -pg
	CPEXT := .DEBUG.LINUX
endif


ifeq ($(TO_BUILD), SHAREDLIB)
	FFLAGS += -fPIC
	LDFLAGS := -shared
	.DEFAULT_GOAL := $(SHAREDLIB)

else ifeq ($(TO_BUILD), STATICLIB)
	.DEFAULT_GOAL := $(STATICLIB)

else ifeq ($(TO_BUILD), EXE)
	# If we are building a standalone executable, input and output can only be via ASCII files
	INPUT_AND_OUTPUT := FILEINOUT
	.DEFAULT_GOAL := $(EXE)
endif


.PHONY := clean
.SECONDARY := $(MOD)
.SECONDARY := $(MAKEFILE)


$(STATICLIB) : $(OBJS)
	$(AR) $(ARFLAGS) $(STATICLIB) $(OBJS)
	cp $(STATICLIB)  ../lib/$(STATICLIB)$(CPEXT).$(INPUT_AND_OUTPUT)
	cp $(HEADER) ../inc/

$(SHAREDLIB) : $(OBJS)
	$(LD) $(LDFLAGS) -o $(SHAREDLIB) $(OBJS)
	cp $(SHAREDLIB)  ../lib/$(SHAREDLIB)$(CPEXT).$(INPUT_AND_OUTPUT)
	cp $(HEADER) ../inc/

$(EXE) : $(OBJS)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJS)

$(MOD) : $(MODSRC)
	$(FC) $(FFLAGS) -D$(INPUT_AND_OUTPUT) -c $(MODSRC)

%.o : %.f03 $(MOD) $(MAKEFILE)
	$(FC) $(FFLAGS) -D$(INPUT_AND_OUTPUT) -D$(TO_BUILD) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(STATICLIB)
	rm -f $(SHAREDLIB)
	rm -f $(EXE)

print-%  : ; @echo $* = $($*)

# To look at symbols in static library when built: nm -C libCShore.a
