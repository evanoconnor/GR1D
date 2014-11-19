include ../../make.inc

NT_SOURCES=nulibtable.F90 \
	   nulibtable_reader.F90 \
	   linterp_many_mod.F90

NT_OBJECTS=$(NT_SOURCES:.F90=.o )

#take care of EOS dependences etc
EXTRAINCS = $(HDF5INCS)
EXTRAOBJECTS = $(HDF5LIBS)

all: nulibtable.a

$(NT_OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

nulibtable.a: $(NT_OBJECTS)
	ar r nulibtable.a *.o

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
