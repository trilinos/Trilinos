include TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include makefile.$(TRILINOS_TARGET)

allbutesipetra: aztecoo epetra ifpack komplex triutils
all: aztecoo epetra epetra_esi ifpack komplex petra triutils

libinstall:
	@echo " "
	@echo " "
	@echo "--------------------------------------------------------------"
	@echo "Building $(CUR_LIB) for the following target:" $(TRILINOS_TARGET)
	@echo "--------------------------------------------------------------"
	@echo " "
	cd $(TRILINOS_HOME)/packages/$(CUR_SRCPATH)$(CUR_LIB)/; mkdir -p obj/$(TRILINOS_TARGET)
	(cd $(TRILINOS_HOME)/packages/$(CUR_SRCPATH)$(CUR_LIB)/obj/$(TRILINOS_TARGET); rm -f Makefile *.[fch]*; $(LN)  ../../src/Makefile . ; $(LN)  ../../src/*.[fch]* .; $(MAKE)) 

anasazi:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
	     libinstall CUR_LIB=anasazi

aztecoo: 
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		libinstall CUR_LIB=aztecoo

epetra:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		libinstall CUR_LIB=epetra

epetra_esi:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
	        libinstall CUR_SRCPATH=../adaptors/ CUR_LIB=epetra_esi

ifpack:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		 libinstall CUR_LIB=ifpack

komplex:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		 libinstall CUR_LIB=komplex

petra:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		libinstall CUR_LIB=petra

tsf:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		 libinstall CUR_LIB=tsf

triutils:
	$(MAKE) TRILINOS_ARCH="$(TRILINOS_ARCH)" TRILINOS_COMM="$(TRILINOS_COMM)"\
		TRILINOS_ID="$(TRILINOS_ID)"\
		 libinstall CUR_LIB=triutils



clean: clean_aztecoo clean_epetra clean_epetra_esi clean_ifpack \
       clean_komplex clean_petra clean_triutils

objclean:
	@echo " "
	@echo " "
	@echo "--------------------------------------------------------------"
	@echo "Removing $(CUR_LIB) object code for the target:" $(TRILINOS_TARGET)
	@echo "--------------------------------------------------------------"
	@echo " "
	rm -rf $(TRILINOS_HOME)/packages/$(CUR_SRCPATH)$(CUR_LIB)/obj/$(TRILINOS_TARGET)

clean_aztecoo:
	$(MAKE) objclean CUR_LIB=aztecoo

clean_epetra:
	$(MAKE) objclean CUR_LIB=epetra

clean_epetra_esi:
	$(MAKE) objclean CUR_SRCPATH=../adaptors/ CUR_LIB=epetra_esi

clean_ifpack:
	$(MAKE) objclean CUR_LIB=ifpack

clean_komplex:
	$(MAKE) objclean CUR_LIB=komplex

clean_petra:
	$(MAKE) objclean CUR_LIB=petra

clean_triutils:
	$(MAKE) objclean CUR_LIB=triutils

clean_tsf:
	$(MAKE) objclean CUR_LIB=tsf

clobber: clobber_aztecoo clobber_epetra clobber_epetra_esi \
         clobber_ifpack clobber_komplex \
         clobber_petra clobber_triutils
	rm -rf $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)

libclobber:
	@echo " "
	@echo " "
	@echo "--------------------------------------------------------------"
	@echo "Removing $(CUR_LIB) library for the target:" $(TRILINOS_TARGET)
	@echo "--------------------------------------------------------------"
	@echo " "
	cd $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET); rm -f lib$(CUR_LIB).a

clobber_aztecoo: clean_aztecoo
	$(MAKE) libclobber CUR_LIB=aztecoo

clobber_epetra: clean_epetra
	$(MAKE) libclobber CUR_LIB=epetra

clobber_epetra_esi: clean_epetra_esi
	$(MAKE) libclobber CUR_LIB=epetra_esi

clobber_ifpack: clean_ifpack
	$(MAKE) libclobber CUR_LIB=ifpack

clobber_komplex: clean_komplex
	$(MAKE) libclobber CUR_LIB=komplex

clobber_petra: clean_petra
	$(MAKE) libclobber CUR_LIB=petra

clobber_triutils: clean_triutils
	$(MAKE) libclobber CUR_LIB=triutils

clobber_tsf: clean_tsf
	$(MAKE) libclobber CUR_LIB=tsf

again: clobber all
