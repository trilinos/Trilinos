HOST=$(shell hostname)
include ../../Makefile.inc

info:
	@echo "make test EXEC=exec or make all"
	@echo " - $(HOST)"
	@echo " - kokkos is located at $(KOKKOS)"
	@echo " - exec is one of $(EXECS)"

EXEC  = KokkosBatched_Test_BlockCrs
EXECS = \
	KokkosBatched_Test_Gemm \
	KokkosBatched_Test_Trsm \
	KokkosBatched_Test_LU \
	KokkosBatched_Test_BlockCrs

# KokkosKernels_Test_Gemm_Complex
# KokkosKernels_Test_Gemv
# 

include $(KOKKOS)/Makefile.kokkos

%.o : %.cpp
	@echo "Compiling $<"
	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) $(INCLUDES) $(KOKKOS_CXXFLAGS) $(KOKKOS_CPPFLAGS) -o $@ -c $< 
#	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) $(INCLUDES) OB$(KOKKOS_CXXFLAGS) $(KOKKOS_CPPFLAGS) -o $@.S -S -c $< 

test : $(EXEC).o
	@echo "Linking $<"
	$(CXX) $(MKL_CFLAG) $(CXXFLAGS) -o $(EXEC) $(EXEC).o $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(LIBS)

all :
	for e in $(EXECS) ; do \
		make test EXEC=$$e; \
	done

clean :
	@/bin/rm -f *.o *~ *.S *.optrpt
	@for e in $(EXECS); do rm -f $$e $$e.o; done
