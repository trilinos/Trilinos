all:
	@echo "Using TRILINOS_HOME=$(TRILINOS_HOME)"
	@echo "Using TRILINOS_ARCH=$(TRILINOS_ARCH)"
	@echo "Using TRILINOS_LIB=$(TRILINOS_LIB)"
	(cd epetra; make; cd -);
	(cd aztecoo; make; cd -);
	(cd ifpack; make; cd -);
	(cd ml; make; cd -);
	(cd nox; make; cd -);
	(cd triutils; make; cd -);
	@echo "NOTE: Amesos examples have not been complied";

clean:
	(cd epetra; make clean; cd -);
	(cd aztecoo; make clean; cd -);
	(cd ifpack; make clean; cd -);
	(cd ml; make clean; cd -);
	(cd nox; make clean; cd -);
	(cd triutils; make clean; cd -);
	(cd triutils; make clean; cd -);
