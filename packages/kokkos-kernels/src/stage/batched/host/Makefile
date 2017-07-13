HOST=$(shell hostname)
include ../Makefile.inc

info:
	@echo "make test EXEC=exec or make all"
	@echo " - $(HOST)"
	@echo " - kokkos is located at $(KOKKOS)"
	@echo " - exec is one of $(EXECS)"

EXEC  = KokkosKernels_Test_BlockCrs
EXECS = \
	KokkosKernels_Test_Gemm \
	KokkosKernels_Test_Trsm \
	KokkosKernels_Test_LU \
	KokkosKernels_Test_BlockCrs

include $(KOKKOS)/Makefile.kokkos

%.o : %.cpp ../*.hpp *.hpp
	@echo "Compiling $<"
	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) $(INCLUDES) $(KOKKOS_CXXFLAGS) $(KOKKOS_CPPFLAGS) -o $@ -c $< 
#	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) $(INCLUDES) OB$(KOKKOS_CXXFLAGS) $(KOKKOS_CPPFLAGS) -o $@.S -S -c $< 

test : $(EXEC).o
	@echo "Linking $<"
	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) -o $(EXEC) $(EXEC).o $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(LIBS)

intel :
	icc $(MKL_CFLAG) gemm_cmpt.c -qopenmp -O3 -o compact_dgemm $(LIBS)
#	icc gemm_cmpt.c -Wl,--start-group ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.a ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.a ${MKL_ROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -I${MKL_ROOT}/include -qopenmp -O3      34 -o compact_dgemm


all :
	for e in $(EXECS) ; do \
		make test EXEC=$$e; \
	done

clean :
	@/bin/rm -f *.o *~ *.S *.optrpt
	@for e in $(EXECS); do rm -f $$e $$e.o; done
