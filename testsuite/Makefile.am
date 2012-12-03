SUBDIRS = test01 test02 test03 test04 test05

TESTS = testset1 testset2 testset3 testset4 testset5 testset6 testset7 testset8 testset9 testset10 testset11

BUILT_SOURCES =  cmpTest
EXTRA_DIST = test.src cmpTest scripts/mktest test1.in test2.in test3.in test4.in test5.in test6.in test7.in test8.in test9.in test10.in test11.in
	 
CLEANFILES = testset1 testset2 testset3 testset4 testset5 testset6 testset7 testset8 testset9 testset10 testset11

################

testset1: setR test1.in
	${TCLSH} cmpTest ${TCLSH} testset1 test1.in ${VERSION}
	
testset2: test2.in
	${TCLSH} cmpTest ${TCLSH} testset2 test2.in ${VERSION}
	
testset3: test3.in
	${TCLSH} cmpTest ${TCLSH} testset3 test3.in ${VERSION}
	
testset4: test4.in
	${TCLSH} cmpTest ${TCLSH} testset4 test4.in ${VERSION}

testset5: test5.in
	${TCLSH} cmpTest ${TCLSH} testset5 test5.in ${VERSION}

testset6: test6.in
	${TCLSH} cmpTest ${TCLSH} testset6 test6.in ${VERSION}

testset7: test7.in
	${TCLSH} cmpTest ${TCLSH} testset7 test7.in ${VERSION}

testset8: test8.in
	${TCLSH} cmpTest ${TCLSH} testset8 test8.in ${VERSION}

testset9: test9.in
	${TCLSH} cmpTest ${TCLSH} testset9 test9.in ${VERSION}
	
testset10: test10.in
	${TCLSH} cmpTest ${TCLSH} testset10 test10.in ${VERSION}
	
testset11: test11.in
	${TCLSH} cmpTest ${TCLSH} testset11 test11.in ${VERSION}
	
#################

setR:
	chmod +x ../src/scripts/nose_${VERSION}
	
