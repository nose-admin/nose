noinst_LIBRARIES = libmolcinp.a libmolc.a
libmolcinp_a_SOURCES = input_molc.F90 resources_molc.F90 
libmolc_a_SOURCES = module_molc.F90 excband.F90
BUILT_SOURCES = input_molc.mod module_molc.mod resources_molc.mod excband.mod

INCLUDES    = -I./ ${LAPACK95_INC} -I${top_srcdir}/lib -I${top_srcdir}/src/pert 
CLEANFILES = ${BUILT_SOURCES}

# suffixes
.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

#
#
#
resources_molc.mod: ../../libresources.a
module_molc.mod: ../../libprepare.a excband.mod
module_molc.o: excband.mod

input_molc.mod: resources_molc.mod 
input_molc.o: resources_molc.mod 


#
#
#
../../libresources.a: 
	cd ../../; make libresources.a
	
../../libprepare.a:
	cd ../../; make libprepare.a

