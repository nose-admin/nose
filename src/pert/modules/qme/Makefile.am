noinst_LIBRARIES = libqmeinp.a libqme.a
libqmeinp_a_SOURCES = input_qme.F90 resources_qme.F90 
libqme_a_SOURCES = module_qme.F90 qme_nakajima_zwanzig.F90 qme_weak_excitons.F90 qme_excitonic_networks.F90 qme_vibronic_networks.F90
BUILT_SOURCES = input_qme.mod module_qme.mod resources_qme.mod qme_nakajima_zwanzig.mod qme_weak_excitons.mod qme_excitonic_networks.mod qme_vibronic_networks.mod

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
resources_qme.mod: ../../libresources.a
qme_nakajima_zwanzig.mod: ../../libprepare.a
qme_weak_excitons.mod: ../../libprepare.a
module_qme.mod: ../../libprepare.a
 
	
module_qme.mod: resources_qme.mod qme_nakajima_zwanzig.mod qme_weak_excitons.mod qme_excitonic_networks.mod qme_vibronic_networks.mod
module_qme.o: resources_qme.mod qme_nakajima_zwanzig.mod qme_weak_excitons.mod qme_excitonic_networks.mod qme_vibronic_networks.mod
input_qme.mod: resources_qme.mod 
input_qme.o: resources_qme.mod 
qme_nakajima_zwanzig.mod: resources_qme.mod 
qme_nakajima_zwanzig.o: resources_qme.mod
qme_weak_excitons.mod: resources_qme.mod 
qme_weak_excitons.o: resources_qme.mod

#
#
#
../../libresources.a: 
	cd ../../; make libresources.a
	
../../libprepare.a:
	cd ../../; make libprepare.a
