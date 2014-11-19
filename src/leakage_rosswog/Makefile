include ../../make.inc

EXTRAINC= -I../nuc_eos -I../

all: leakage_rosswog.a

leakage_rosswog.a: leak_rosswog.o
	ar r leakage_rosswog.a leak_rosswog.o

leak_rosswog.o: leak_rosswog.F90
	$(F90) $(F90FLAGS) -c leak_rosswog.F90 $(EXTRAINC)

clean:
	rm -f *.o *.mod *.a
