DIRS = CCA_Wrapper TestComponents

default : component test

component :
	echo " Making in CCA_Wrapper ... " ; \
	cd CCA_Wrapper/ ; make ; cd .. 

static_lib :
	echo " Making in  CCA_Wrapper ... " ; \
	cd CCA_Wrapper ; make -f Makefile.static ; cd ..

test :
	echo " Making in TestComponents ... " ; \
	cd TestComponents/; make ; cd .. 

test_static_lib :
	echo " Making in TestComponents ... " ; \
	cd TestComponents/ make -f Makefile.static ; cd .. 

clean :
	for dir in $(DIRS) ; \
	do echo " Making clean in $$dir ..." ;\
	cd $$dir ; make -f Makefile.static clean ; make clean ; cd .. ;\
	done

distclean :
	for dir in $(DIRS) ; \
	do echo " Making distclean in $$dir ..." ;\
	cd $$dir ; make distclean ; cd .. ;\
	done
	/bin/rm -fr config.log config.status autom4te* MakeIncl.ZoltanComponent
