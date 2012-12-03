noinst_LIBRARIES = libtestinp.a libtest.a
libtestinp_a_SOURCES = input_test.F90 resources_test.F90
libtest_a_SOURCES = module_test.F90 
BUILT_SOURCES = input_test.mod module_test.mod resources_test.mod 

INCLUDES    = -I./ -I${top_srcdir}/lib -I${top_srcdir}/src/pert
CLEANFILES = ${BUILT_SOURCES}

# suffixes
.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

#
#
#
resources_test.mod: ../../libresources.a
module_test.mod: ../../libprepare.a

input_test.mod: resources_test.mod 
input_test.o: resources_test.mod 

#
#
#
../../libresources.a: 
	cd ../../; make libresources.a
	
../../libprepare.a:
	cd ../../; make libprepare.a

