#=======================================================================
# create ML library file and executable
#=======================================================================

DEFAULT:
	echo "Making ML library"
	@(cp -f Coarsen/*.h Include)
	@(cp -f Comm/*.h Include)
	@(cp -f FEGrid/*.h Include)
	@(cp -f Krylov/*.h Include)
	@(cp -f Main/*.h Include)
	@(cp -f Operator/*.h Include)
	@(cp -f Smoother/*.h Include)
	@(cp -f Utils/*.h Include)
	@(cd Include; make)
	@(cd Coarsen; make)
	@(cd Comm; make)
	@(cd FEGrid; make)
	@(cd Krylov; make)
	@(cd Main; make)
	@(cd Operator; make)
	@(cd Smoother; make)
	@(cd Utils; make)
	@(cd lib; make)

clean:
	@(cd Coarsen; make clean)
	@(cd Comm; make clean)
	@(cd FEGrid; make clean)
	@(cd Krylov; make clean)
	@(cd Main; make clean)
	@(cd Operator; make clean)
	@(cd Smoother; make clean)
	@(cd Utils; make clean)
	@(cd lib; make clean)
	@(cd Include; make clean)
	@(cd Runtime; make clean)

