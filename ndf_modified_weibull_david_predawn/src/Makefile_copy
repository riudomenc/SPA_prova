###### NETCDF SETUP #######
#
NC_LIB_DIR=/usr/lib
NC_INC_DIR=/usr/include
LAPACK_LIB_DIR=/usr/lib/lapack
BLAS_LIB_DIR=/usr/lib/libblas
#
NCFLAGS=-I$(NC_INC_DIR) -L$(NC_LIB_DIR) -L$(LAPACK_LIB_DIR) -L$(BLAS_LIB_DIR) -lnetcdff -lnetcdf -llapack -lblas

###############################

###### Compiler choice #####
#
# to use the gnu compiler..
#
FC = gfortran
#for normal...
#FCFLAGS = -cpp -ffree-line-length-none $(NCFLAGS)
#for debugging...
FCFLAGS = -cpp -fbounds-check -fimplicit-none -ffree-line-length-none -frange-check -ftree-vectorizer-verbose=0 -ggdb -O0 -Wall $(NCFLAGS) -std=gnu -pedantic
#for profiling...FCFLAGS = -cpp -p -g -ffree-line-length-none $(NCFLAGS)
#
# to use the intel compiler...
#
#FC = ifort
#for normal......
#FCFLAGS = -fpp -vec-report=0 $(NCFLAGS)
#for debugging...FCFLAGS = -fpp -check all -debug-parameters all -ftrapuv -g -fpe0 -implicitnone -O0 -p -stand f03 -traceback -vec-report=0 $(NCFLAGS) 
#for profiling...FCFLAGS = -fpp -p -g $(NCFLAGS)
#
##############################

spa_SOURCES = enkf.f90  enkf_declarations.f90  all_declarations.f90  allocate_carbon.f90  canopy.f90  config_tools.f90 \
              leaf.f90  light.f90  log_tools.f90  main.f90  math_tools.f90  netcdf_tools.f90  \
              scale_declarations.f90  soil_air.f90  soil_functions.f90  spa_cmd_line.f90 spa_config.f90  \
              spa_initialise.f90  spa_io.f90  spa_io_csv.f90  spa_io_netcdf.f90

spa:	$(spa_SOURCES:%.f90=%.o)
	$(FC) -i-dynamic $(NCFLAGS) -o $@ $^
	mv spa ../

clean:
	rm -f *.mod *.o ../spa make_dep

%.o:	%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

Makefile:	make_dep

make_dep:	$(spa_SOURCES)
		./f90_dependencies.py $^ > $@

-include make_dep
