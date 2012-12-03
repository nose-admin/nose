SUBDIRS = modules

# executable file destination
binnosedir = ${bindir}/nose-home/${VERSION}/bin

##
# NOSE Driver
##
binnose_PROGRAMS = ndrv.x 
ndrv_x_SOURCES =   main.F90 modules.F90 main_comm.F90 

ndrv_x_LDADD = -L./ -L${top_srcdir}/lib/ \
 -L${top_srcdir}/src/pert/modules/tdpt3 \
 -L${top_srcdir}/src/pert/modules/qme \
 -L${top_srcdir}/src/pert/modules/montecarlo \
 -L${top_srcdir}/src/pert/modules/molc \
 -L${top_srcdir}/src/pert/modules/test \
 -L${top_srcdir}/src/pert \
 -ltdpt3 -lqme -lmontecarlo -lmolc -ltest \
 -lprepare \
 -ltdpt3inp -lqmeinp -lmontecarloinp -lmolcinp -ltestinp \
 -lresources \
 -lnose  ${LAPACK_LIBS} ${BLAS_LIBS} ${MPILIBS} # -lgfortran


#
# dependencies on NOSE module files
#
ndrv_x_DEPENDENCIES = modules/tdpt3/libtdpt3.a modules/qme/libqme.a modules/test/libtest.a modules/montecarlo/libmontecarlo.a modules/molc/libmolc.a

noinst_LIBRARIES = libresources.a libprepare.a

libresources_a_SOURCES = resources.F90 nis.F90
libprepare_a_SOURCES = prepare.F90 bath_interactions.F90 excitons.F90 \
 simulate_wizard.F90 goft.F90 rates.F90 input_wizard.F90 orfact.F90 resp3_rss.F90 resp3_wd.F90 resp3_site.F90 resp3_ext.F90 \
 response3.F90 signal3.F90 twod.F90 pump_probe.F90 nakajima_zwanzig_shared.F90
 
 
#
# dependencies of NOSE module files
#
libprepare_a_DEPENDENCIES = modules/tdpt3/libtdpt3inp.a modules/qme/libqmeinp.a modules/test/libtestinp.a modules/montecarlo/libmontecarloinp.a modules/molc/libmolcinp.a


INCLUDES = -I. -I${top_srcdir}/lib/ -I${top_srcdir}/src/pert/modules/tdpt3 \
  -I${top_srcdir}/src/pert/modules/qme -I${top_srcdir}/src/pert/modules/test -I${top_srcdir}/src/pert \
  -I${top_srcdir}/src/pert/modules/montecarlo \
  -I${top_srcdir}/src/pert/modules/molc
  
BUILT_SOURCES = resources.mod nis.mod libresources.a modules/tdpt3/libtdpt3inp.a modules/qme/libqmeinp.a \
  modules/montecarlo/libmontecarloinp.a \
  modules/test/libtestinp.a \
  bath_interactions.mod \
  excitons.mod prepare.mod  input_wizard.mod modules.mod \
  goft.mod rates.mod  simulate_wizard.mod orfact.mod resp3_rss.mod response3.mod resp3_wd.mod resp3_site.mod resp3_ext.mod \
  signal3.mod twod.mod pump_probe.mod  nakajima_zwanzig_shared.mod \
  modules/tdpt3/libtdpt3.a \
  modules/qme/libqme.a modules/test/libtest.a modules/montecarlo/libmontecarlo.a modules/molc/libmolc.a\
  main_comm.mod

CLEANFILES = ${BUILT_SOURCES}

LNDEP = ${top_srcdir}/lib/libnose.a

FC = ${MPIFC}

# suffixes
.F90.mod:
	${MPIFC} ${FCFLAGS} ${INCLUDES}  -c $<

.F90.o:
	${MPIFC} ${FCFLAGS} ${INCLUDES}  -c $<


	
##################################################################################
#                      NOSE LIBRARY and EXTERNAL LIBRARIES
##################################################################################

${top_srcdir}/lib/libnose.a:
	cd ${top_srcdir}/lib; make all
	
${top_srcdir}/lib/std_types.mod:
	cd ${top_srcdir}/lib; make std_types.mod

##################################################################################
#                            END NOSE LIBRARY
##################################################################################
#
##################################################################################
#                            MODULE LIBRARIES
##################################################################################

##
#  QME
##

modules/qme/libqmeinp.a: libresources.a
	cd modules/qme; make libqmeinp.a

modules/qme/libqme.a: libprepare.a
	cd modules/qme; make libqme.a


##
# TDPT3
##

modules/tdpt3/libtdpt3inp.a: libresources.a
	cd modules/tdpt3; make libtdpt3inp.a

modules/tdpt3/libtdpt3.a: libprepare.a
	cd modules/tdpt3; make libtdpt3.a

##
# MOLC
##

modules/molc/libmolcinp.a: libresources.a
	cd modules/molc; make libmolcinp.a

modules/molc/libmolc.a: libprepare.a
	cd modules/molc; make libmolc.a


##
# TEST
##

modules/test/libtestinp.a: libresources.a
	cd modules/test; make libtestinp.a

modules/test/libtest.a: libprepare.a
	cd modules/test; make libtest.a


##
#  MC
##

modules/montecarlo/libmontecarloinp.a: libresources.a
	cd modules/montecarlo; make libmontecarloinp.a

modules/montecarlo/libmontecarlo.a: libprepare.a
	cd modules/montecarlo; make libmontecarlo.a

##################################################################################
#                           END MODULE LIBRARIES
##################################################################################
#
##################################################################################
#                          NOSE DRIVER DEPENDENCIES
##################################################################################

nis.mod: ${top_srcdir}/lib/std_types.mod ${LNDEP}
resources.mod: nis.mod 

##
#   connection point 1 for the modules
##
simulate_wizard.mod: resources.mod modules/tdpt3/libtdpt3inp.a modules/qme/libqmeinp.a modules/test/libtestinp.a modules/montecarlo/libmontecarloinp.a


##
#   connection point 2 for the modules
## 

input_wizard.mod:  simulate_wizard.mod   
excitons.mod: input_wizard.mod 
goft.mod: excitons.mod 
rates.mod: goft.mod  
bath_interactions.mod: goft.mod rates.mod
prepare.mod: bath_interactions.mod orfact.mod
orfact.mod: resources.mod
resp3_rss.mod:
resp3_wd.mod: resp3_rss.mod 
resp3_site.mod: resp3_rss.mod 
resp3_ext.mod: resp3_rss.mod 
response3.mod: resp3_rss.mod resp3_wd.mod resp3_site.mod resp3_ext.mod resources.mod
signal3.mod: response3.mod
twod.mod: signal3.mod
nakajima_zwanzig_shared.mod: 

##
#    connection point 3 for the modules
##
modules.mod: modules/tdpt3/libtdpt3.a modules/qme/libqme.a modules/test/libtest.a modules/montecarlo/libmontecarlo.a modules/molc/libmolc.a
main_comm.mod: prepare.mod modules.mod
main.o: main_comm.mod

##################################################################################
#                           END NOSE DRIVER DEPENDECIES
##################################################################################