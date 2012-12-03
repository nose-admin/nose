noinst_LIBRARIES = libtdpt3inp.a libtdpt3.a
libtdpt3inp_a_SOURCES = input_tdpt3.F90 resources_tdpt3.F90
libtdpt3_a_SOURCES = module_tdpt3.F90 response.F90
BUILT_SOURCES = input_tdpt3.mod module_tdpt3.mod resources_tdpt3.mod response.mod

INCLUDES    = -I./ ${LAPACK95_INC} -I${top_srcdir}/lib -I${top_srcdir}/src/pert -I${top_srcdir}/src/pert/modules/qme
CLEANFILES = ${BUILT_SOURCES}

# suffixes
.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

#
#
#
resources_tdpt3.mod: ../../libresources.a
module_tdpt3.mod: ../../libprepare.a

input_tdpt3.mod: resources_tdpt3.mod 
input_tdpt3.o: resources_tdpt3.mod 

#
#
#
../../libresources.a: 
	cd ../../; make libresources.a
	
../../libprepare.a:
	cd ../../; make libprepare.a

