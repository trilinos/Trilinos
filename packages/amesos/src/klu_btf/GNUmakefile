# for testing only:
# TEST = -DTESTING

INCLUDES =  -I../klu_kernel -I../AMD/Include -I../btf -I../harwell -I../AMD/Source -I../colamd
CMEX = mex -inline $(INCLUDES) $(TEST)
FMEX = mex $(INCLUDES)

all: klua.mexglx kluf.mexglx klus.mexglx klua_harwell.mexglx \
    klus_harwell.mexglx kluag.mexglx klust.mexglx

AMD = amd_aat.o amd_1.o amd_2.o amd_dump.o amd_postorder.o amd_post_tree.o \
    amd_defaults.o amd_order.o amd_control.o amd_info.o amd_valid.o

COLAMD = colamd.o
COLAMDH = ../colamd/colamd.h

KLU_BTF = klu_dump.o klu_btf_free_symbolic.o \
    klu_btf_factor.o klu_btf_free_numeric.o klu_btf_solve.o klu_btf_scale.o \
    klu_btf_refactor.o klu_btf_defaults.o klu_btf_analyze_given.o \
    klu_btf_tsolve.o

BTF = btf_order.o maxtrans.o strongcomp.o

KLU = klu.o klu_kernel.o

FHARWELL = harwell.o mc13d.o mc13e.o mc21a.o mc21b.o

CHARWELL = charwell.o

OBJ = $(AMD) $(KLU_BTF) $(KLU) $(FHARWELL) $(CHARWELL) $(COLAMD) $(BTF)

AMDH = ../AMD/Source/amd_internal.h ../AMD/Include/amd.h

INC = klu_btf.h $(AMDH) ../harwell/charwell.h ../klu_kernel/klu.h \
    ../klu_kernel/klu_kernel.h ../klu_kernel/tictoc.h klu_dump.h $(COLAMDH) \
    klu_btf_internal.h ../btf/btf.h ../btf/btf_internal.h \
    ../StandAlone/dsecnd.h

amd_%.o: ../AMD/Source/amd_%.c $(AMDH)
	$(CMEX) -DDINT -c $<

colamd.o: ../colamd/colamd.c $(COLAMDH)
	$(CMEX) -c $<

klu_dump.o: klu_dump.c $(INC)
	$(CMEX) -c $<

klu_btf%.o: klu_btf%.c $(INC)
	$(CMEX) -c $<

klu.o: ../klu_kernel/klu.c $(INC)
	$(CMEX) -c $<

klu_kernel.o: ../klu_kernel/klu_kernel.c $(INC)
	$(CMEX) -c $<

tictoc.o: ../klu_kernel/tictoc.c $(INC)
	$(CMEX) -c $<

mc%.o: ../harwell/mc%.f
	$(FMEX) -c $<

harwell.o: ../harwell/harwell.f
	$(FMEX) -c $<

charwell.o: ../harwell/charwell.c $(INC)
	$(CMEX) -c $<

klua.mexglx: klua_mex.c $(OBJ) $(INC) klu_btf_analyze.c
	$(CMEX) -output klua klua_mex.c $(OBJ) klu_btf_analyze.c

kluag.mexglx: kluag_mex.c $(OBJ) $(INC) klu_btf_analyze.c
	$(CMEX) -output kluag kluag_mex.c $(OBJ) klu_btf_analyze.c

klua_harwell.mexglx: klua_mex.c $(OBJ) $(INC) klu_btf_analyze.c
	$(CMEX) -DHARWELL -output klua_harwell klua_mex.c $(OBJ) klu_btf_analyze.c

kluf.mexglx: kluf_mex.c $(OBJ) $(INC)
	$(CMEX) -output kluf kluf_mex.c $(OBJ)

klus.mexglx: klus_mex.c $(OBJ) $(INC) klu_btf_analyze.c ../StandAlone/dsecnd.c
	$(CMEX) -output klus klus_mex.c $(OBJ) klu_btf_analyze.c ../StandAlone/dsecnd.c

klust.mexglx: klust_mex.c $(OBJ) $(INC) klu_btf_analyze.c ../StandAlone/dsecnd.c
	$(CMEX) -output klust klust_mex.c $(OBJ) klu_btf_analyze.c ../StandAlone/dsecnd.c

klus_harwell.mexglx: klus_mex.c $(OBJ) $(INC) klu_btf_analyze.c ../StandAlone/dsecnd.c
	$(CMEX) -DHARWELL -output klus_harwell klus_mex.c $(OBJ) klu_btf_analyze.c ../StandAlone/dsecnd.c

btf_order.o: ../btf/btf_order.c
	$(CMEX) -c $<

maxtrans.o: ../btf/maxtrans.c
	$(CMEX) -c $<

strongcomp.o: ../btf/strongcomp.c
	$(CMEX) -c $<

clean:
	- rm *.o

purge: clean
	- rm *.mex* *.dll

