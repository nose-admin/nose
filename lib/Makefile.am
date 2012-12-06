SUBDIRS = cuda

noinst_LIBRARIES = libnose.a 

libnose_a_SOURCES = std_types.F90 std_lapack.F90 numer_fft.F90 numer_fft_row.F90 util_char_graphics.F90 util_space.F90 \
                    std_io.F90 numer_ode.F90 numer_ode_rk4.F90 data_list.F90 util_timing.F90 numer_matrix.F90 \
                    numer_interp.F90 sci_foerster.F90 sci_brownian_oscillator.F90 numer_random.F90 std_parallel.F90 sci_redfield.F90 \
                    sci_liouville.F90 sci_misc.F90 util_allocation.F90 util_allocation.h qch_lib.F90 util_tracing.h numer_gamma.F90

BUILT_SOURCES = std_types.mod std_lapack.mod numer_fft.mod util_char_graphics.mod \
 util_space.mod std_io.mod numer_ode.mod data_list.mod sci_misc.mod \
 util_timing.mod numer_matrix.mod numer_interp.mod sci_foerster.mod sci_brownian_oscillator.mod \
 numer_random.mod  std_parallel.mod sci_redfield.mod sci_liouville.mod util_allocation.mod qch_lib.mod \
 numer_gamma.mod 


CLEANFILES = ${BUILT_SOURCES}  


ALLDEP = std_types.mod util_allocation.mod

FC = ${MPIFC}

# suffixes
.f90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} -I. -c $<

.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} -I. -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES} -I. -c $<
	


#
# Dependencies
#

# numer
#+aurelia
numer_gamma.mod: ${ALLDEP} 
numer_gamma.o: ${ALLDEP}
#-aurelia
numer_fft.mod: ${ALLDEP} std_io.mod
numer_fft.o:   ${ALLDEP} std_io.mod
numer_fft_row.o:   ${ALLDEP} std_io.mod
numer_random.mod: ${ALLDEP} numer_interp.mod
numer_random.o: ${ALLDEP} numer_interp.mod
numer_ode.mod: ${ALLDEP} std_io.mod
numer_ode.o: ${ALLDEP} std_io.mod
numer_ode_rk4.o: ${ALLDEP} std_io.mod
numer_matrix.mod: ${ALLDEP} std_io.mod
numer_matrix.o: ${ALLDEP} std_io.mod
numer_interp.mod: ${ALLDEP} numer_matrix.mod
numer_interp.o: ${ALLDEP} numer_matrix.mod

# std
std_io.mod: ${ALLDEP}
std_io.o: ${ALLDEP}
std_lapack.mod: ${ALLDEP}
std_lapack.o: ${ALLDEP}

# util
util_char_graphics.mod: ${ALLDEP} std_io.mod
util_char_graphics.o:   ${ALLDEP} std_io.mod

util_space.mod: ${ALLDEP}
util_space.o:   ${ALLDEP}

util_timing.mod: ${ALLDEP}
util_timing.o: ${ALLDEP}

util_allocation.mod: std_types.mod
util_allocation.o: std_types.mod

# data
data_list.mod: ${ALLDEP}
data_list.o: ${ALLDEP}

# science

sci_foerster.mod: ${ALLDEP} numer_interp.mod
sci_foerster.o: ${ALLDEP} numer_interp.mod

sci_redfield.mod: ${ALLDEP} numer_interp.mod
sci_redfield.o: ${ALLDEP} numer_interp.mod

sci_brownian_oscillator.mod: ${ALLDEP} numer_interp.mod
sci_brownian_oscillator.o: ${ALLDEP} numer_interp.mod

sci_liouville.mod: ${ALLDEP} std_io.mod
sci_liouville.o: ${ALLDEP} std_io.mod

sci_misc.mod: ${ALLDEP} std_io.mod
sci_misc.o: ${ALLDEP} std_io.mod

std_parallel.o: ${ALLDEP} std_io.mod 
std_parallel.mod: ${ALLDEP} std_io.mod 

qch_lib.o: ${ALLDEP} std_types.mod std_io.mod numer_matrix.mod sci_misc.mod
qch_lib.mod: ${ALLDEP} std_types.mod std_io.mod numer_matrix.mod sci_misc.mod





