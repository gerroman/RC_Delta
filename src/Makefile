GCC=clang

DEBUGFLAG=-O0 -g -Wall -pedantic -pedantic-errors -W
RELEASEFLAG=-O3 -Wall -Werror -pedantic -pedantic-errors -W

ROOTCFLAGS=`root-config --cflags`
ROOTLDCFLAGS=`root-config --ldflags`
ROOTLIBS=`root-config --libs`

CFLAGS=$(RELEASEFLAG) $(ROOTCFLAGS)
LFLAGS=-lm -lstdc++

MATRIXOBJS=iMdummy2.o iMe2.o iMp2.o iMpiMe.o iMd2.o iMdiMe.o iMdiMeSoft.o iMdiMeHard.o iMdiMp.o iMborn2.o delta_formfactors.o proton_formfactors.o
INTEGRATEOBJS=integrate.o integrand.o born.o kinematic.o

.PHONY=clean all

all: main.out test_matrix_elements.out

main.out: main.o $(INTEGRATEOBJS) $(MATRIXOBJS)
	$(GCC) $^ $(LFLAGS) $(ROOTLDCFLAGS) $(ROOTLIBS) -o $@

test_matrix_elements.out: test_matrix_elements.o $(MATRIXOBJS)
	$(GCC) $^ $(LFLAGS) -o $@

.cpp.o:
	$(GCC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o *.out
