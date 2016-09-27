FORTRAN = gfortran
SOURCES = driver.f UMAT.f DSID_EISPACK.f DSID_INVERSE.f
DFLAG = -g
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = mTest.x

all : $(EXECUTABLE) $(OBJECTS)

debug :
	${FORTRAN} ${DFLAG} ${SOURCES} -o ${EXECUTABLE}

$(EXECUTABLE) : $(OBJECTS)
	$(FORTRAN) $(OBJECTS) -o $@ 

$(OBJECTS):
	$(FORTRAN) $(SOURCES) -c   

clean :
	rm -f *.out $(EXECUTABLE) *.o

