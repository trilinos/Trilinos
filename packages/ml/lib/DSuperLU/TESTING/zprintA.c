#include "zsp_defs.h"

int zprintA(int n, int lda, doublecomplex *a)
{
    int i, j;
    printf("ar=[\n");
    for (i = 0; i < n; ++i) {
	for (j = 0; j < n; ++j) printf("%16.8e", a[i+j*lda].r);
	printf("\n");
    }
    printf("];\n");
    printf("ai=[\n");
    for (i = 0; i < n; ++i) {
	for (j = 0; j < n; ++j) printf("%16.8e", a[i+j*lda].i);
	printf("\n");
    }
    printf("];\n");
}
