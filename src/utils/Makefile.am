# executable file destination
binnosedir = ${bindir}/nose-home/${VERSION}/bin
##
# NOSE Driver
##
binnose_PROGRAMS = dip.x
dip_x_SOURCES = dipolyf.F90

dip_x_LDADD = -L./ -L${top_srcdir}/lib/ -L${top_srcdir}/src/pert/ \
 -lnose  ${LAPACK_LIBS} ${BLAS_LIBS} ${MPILIBS} # -lgfortran


INCLUDES = -I. -I${top_srcdir}/lib/ -I${top_srcdir}/src/pert/modules

BUILT_SOURCES = 

CLEANFILES = ${BUILT_SOURCES}

LNDEP = ${top_srcdir}/lib/libnose.a

FC = ${MPIFC}

# suffixes
.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES}  -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES}  -c $<
	
