

binnosedir = ${bindir}/nose-home/${VERSION}/bin
datadir    = ${bindir}/nose-home/scilab
confdir    = ${bindir}/nose-home/
conflocdir = ${bindir}/nose-home/

conf_SCRIPTS = nose_${VERSION}
binnose_SCRIPTS = nose2
EXTRA_DIST = nose.tcl ubernose.tcl cmpNose


BUILT_SOURCES = nose.tcl ubernose.tcl cmpNose
	 
CLEANFILES = nose_${VERSION} nose2

#
# "compile" nose script
#
nose_${VERSION}: ubernose.tcl
	${TCLSH} cmpNose ${TCLSH} ${USEMPI} nose_${VERSION} ubernose.tcl ${confdir} ${binnosedir} "${FCFLAGS}"

nose2: nose.tcl
	${TCLSH} cmpNose ${TCLSH} ${USEMPI} nose2 nose.tcl ${confdir} ${binnosedir} "${FCFLAGS}"
	
	

 