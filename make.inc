F90=gfortran
F90FLAGS= -g -O3 -fopenmp
LDFLAGS= -g -O3 -fopenmp 
MODINC="-I ./"

HDF5DIR=/home/x_aurbe/hdf5-1.8.20-mpif90
HDF5INCS=-I${HDF5DIR}/include
HDF5LIBS=-L${HDF5DIR}/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz -DH5_USE_16_API

LAPACKDIR=/proj/astro_extreme/evanoc/opt/lapack-3.8.0
LAPACKLIBS=-L${LAPACKDIR} -llapack -lblas


HAVE_NUC_EOS=1
HAVE_RESTART=1
HAVE_LEAK_ROS=0
HAVE_LAPACK=1
