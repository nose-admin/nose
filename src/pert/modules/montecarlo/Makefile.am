noinst_LIBRARIES = libmontecarloinp.a libmontecarlo.a
libmontecarloinp_a_SOURCES = input_montecarlo.F90 resources_montecarlo.F90
libmontecarlo_a_SOURCES = module_montecarlo.F90
BUILT_SOURCES = input_montecarlo.mod module_montecarlo.mod resources_montecarlo.mod

INCLUDES    = -I./ ${LAPACK95_INC} -I${top_srcdir}/lib -I${top_srcdir}/src/pert -I${top_srcdir}/src/pert/modules/tdpt3
CLEANFILES = ${BUILT_SOURCES}

# suffixes

.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES} ${LAPACK95_INC} -c $<
#
#
#
resources_montecarlo.mod: ../../libresources.a
module_montecarlo.mod: ../../libresources.a
module_montecarlo.mod: ../../libprepare.a
 
	
module_montecarlo.mod: resources_montecarlo.mod
module_montecarlo.o: resources_montecarlo.mod module_montecarlo.mod
#module_qme.mod: resources_qme.mod qme_nakajima_zwanzig.mod qme_weak_excitons.mod
#module_qme.o: resources_qme.mod qme_nakajima_zwanzig.mod qme_weak_excitons.mod
input_montecarlo.mod: resources_montecarlo.mod 
input_montecarlo.o: resources_montecarlo.mod 
#qme_nakajima_zwanzig.mod: resources_qme.mod 
#qme_nakajima_zwanzig.o: resources_qme.mod
#qme_weak_excitons.mod: resources_qme.mod 
#qme_weak_excitons.o: resources_qme.mod

#
#
#
../../libresources.a: 
	cd ../../; make libresources.a
	
../../libprepare.a:
	cd ../../; make libprepare.a
