
all: maxtrans.mexglx maxtrans_recursive.mexglx strongcomp.mexglx \
	strongcomp_recursive.mexglx btf.mexglx

just: strongcomp.mexglx strongcomp_recursive.mexglx

maxtrans.mexglx: maxtrans.c btf.h maxtrans_mex.c btf_internal.h
	mex maxtrans.c maxtrans_mex.c

maxtrans_recursive.mexglx: maxtrans.c btf.h maxtrans_mex.c btf_internal.h
	mex -DRECURSIVE -output maxtrans_recursive maxtrans.c maxtrans_mex.c

strongcomp.mexglx: strongcomp.c btf.h strongcomp_mex.c btf_internal.h
	mex strongcomp.c strongcomp_mex.c

strongcomp_recursive.mexglx: strongcomp.c btf.h strongcomp_mex.c btf_internal.h
	mex -DRECURSIVE -output strongcomp_recursive strongcomp.c strongcomp_mex.c


btf.mexglx: strongcomp.c btf.h btf_internal.h maxtrans.c btf_mex.c btf_order.c
	mex -output btf maxtrans.c strongcomp.c btf_mex.c btf_order.c

purge:
	rm *.o *.mex*
