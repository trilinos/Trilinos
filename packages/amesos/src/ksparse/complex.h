/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 * Copyright (c) 1985 Thomas L. Quarles
 */
#ifndef CMPLX
#define CMPLX "complex.h $Revision$  on $Date$ "

/*  header file containing definitions for complex functions
 *
 *  Each expects two arguments for each complex number - a real and an
 *  imaginary part.
 */
typedef struct {
    double real;
    double imag;
} SPcomplex;


#define DC_ABS(a,b) (FABS(a) + FABS(b))

#ifdef notdef
#define DC_DIV(a,b,c,d,x,y) { \
    double r,s;\
    if(FABS(c)>FABS(d)) { \
        r=(d)/(c);\
        s=(c)+r*(d);\
        x=((a)+(b)*r)/s;\
        y=((b)-(a)*r)/s;\
    } else { \
        r=(c)/(d);\
        s=(d)+r*(c);\
        x=((a)*r+(b))/s;\
        y=((b)*r-(a))/s;\
    }\
}
#endif /*notdef */

#ifndef HAS_SHORTMACRO
#define DC_DIVEQ(a,b,c,d) { \
    double r,s,x,y;\
    if(FABS(c)>FABS(d)) { \
        r=(d)/(c);\
        s=(c)+r*(d);\
        x=((*(a))+(*(b))*r)/s;\
        y=((*(b))-(*(a))*r)/s;\
    } else { \
        r=(c)/(d);\
        s=(d)+r*(c);\
        x=((*(a))*r+(*(b)))/s;\
        y=((*(b))*r-(*(a)))/s;\
    }\
    (*(a)) = x; \
    (*(b)) = y; \
}
#else /*HAS_SHORTMACRO*/
#define DC_DIVEQ DCdiveq
#ifdef __STDC__
extern void DCdiveq(double*,double*,double,double);
#else /* stdc */
extern void DCdiveq();
#endif /* stdc */
#endif /*HAS_SHORTMACRO*/

#ifndef HAS_SHORTMACRO
#define DC_MULT(a,b,c,d,x,y) { \
    *(x) = (a) * (c) - (b) * (d) ;\
    *(y) = (a) * (d) + (b) * (c) ;\
}
#else /*HAS_SHORTMACRO*/
#define DC_MULT DCmult
#ifdef __STDC__
extern void DCmult(double,double,double,double,double*,double*);
#else /* stdc */
extern void DCmult();
#endif /* stdc */
#endif /*HAS_SHORTMACRO*/

#ifdef notdef
#define DC_MINUS(a,b,c,d,x,y) { \
    (x) = (a) - (c) ;\
    (y) = (b) - (d) ;\
}
#endif /*notdef*/

#ifndef HAS_SHORTMACRO
#define DC_MINUSEQ(a,b,c,d) { \
    *(a) -= (c) ;\
    *(b) -= (d) ;\
}
#else /*HAS_SHORTMACRO*/
#define DC_MINUSEQ DCminusEq
#ifdef __STDC__
extern void DCminusEq(double*,double*,double,double);
#else /* stdc */
extern void DCminusEq();
#endif /* stdc */
#endif /*HAS_SHORTMACRO*/

#define	C_SQRT(A) {							      \
	double	_mag, _a;						      \
	if ((A).imag == 0.0) {						      \
	    if ((A).real < 0.0) {					      \
		(A).imag = sqrt(-(A).real);				      \
		(A).real = 0.0;						      \
	    } else {							      \
		(A).real = sqrt((A).real);				      \
		(A).imag = 0.0;						      \
	    }								      \
	} else {							      \
	    _mag = sqrt((A).real * (A).real + (A).imag * (A).imag);	      \
	    _a = (_mag - (A).real) / 2.0;				      \
	    if (_a <= 0.0) {						      \
		(A).real = sqrt(_mag);					      \
		(A).imag /= (2.0 * (A).real); /*XXX*/			      \
	    } else {							      \
		_a = sqrt(_a);						      \
		(A).real = (A).imag / (2.0 * _a);			      \
		(A).imag = _a;						      \
	    }								      \
	}								      \
    }

#define	C_MAG2(A) (((A).real = (A).real * (A).real + (A).imag * (A).imag),    \
	(A).imag = 0.0)

#define	C_CONJ(A) ((A).imag *= -1.0)

#define	C_CONJEQ(A,B) {							      \
	(A).real = (B.real);						      \
	(A).imag = - (B.imag);						      \
    }

#define	C_EQ(A,B) {							      \
	(A).real = (B.real);						      \
	(A).imag = (B.imag);						      \
    }

#define	C_NORM(A,B) {							      \
	if ((A).real == 0.0 && (A).imag == 0.0) {			      \
	    (B) = 0;							      \
	} else {							      \
	    while (FABS((A).real) > 1.0 || FABS((A).imag) > 1.0) {	      \
		(B) += 1;						      \
		(A).real /= 2.0;					      \
		(A).imag /= 2.0;					      \
	    }								      \
	    while (FABS((A).real) <= 0.5 && FABS((A).imag) <= 0.5) {	      \
		(B) -= 1;						      \
		(A).real *= 2.0;					      \
		(A).imag *= 2.0;					      \
	    }								      \
	}								      \
    }

#define	C_ABS(A) (sqrt((A).real * (A.real) + (A.imag * A.imag)))

#define	C_MUL(A,B) {							      \
	double	TMP1, TMP2;						      \
	TMP1 = (A.real);						      \
	TMP2 = (B.real);						      \
	(A).real = TMP1 * TMP2 - (A.imag) * (B.imag);			      \
	(A).imag = TMP1 * (B.imag) + (A.imag) * TMP2;			      \
    }

#define	C_MULEQ(A,B,C) {						      \
	(A).real = (B.real) * (C.real) - (B.imag) * (C.imag);		      \
	(A).imag = (B.real) * (C.imag) + (B.imag) * (C.real);		      \
    }

#define	C_DIV(A,B) {							      \
	double	_tmp, _mag;						      \
	_tmp = (A.real);						      \
	(A).real = _tmp * (B.real) + (A).imag * (B.imag);		      \
	(A).imag = - _tmp * (B.imag) + (A.imag) * (B.real);		      \
	_mag = (B.real) * (B.real) + (B.imag) * (B.imag);		      \
	(A).real /= _mag;						      \
	(A).imag /= _mag;						      \
    }

#define	C_DIVEQ(A,B,C) {						      \
	double	_mag;							      \
	(A).real = (B.real) * (C.real) + (B.imag) * (C.imag);		      \
	(A).imag = (B.imag) * (C.real) - (B.real) * (C.imag) ;		      \
	_mag = (C.real) * (C.real) + (C.imag) * (C.imag);		      \
	(A).real /= _mag;						      \
	(A).imag /= _mag;						      \
    }

#define	C_ADD(A,B) {							      \
	(A).real += (B.real);						      \
	(A).imag += (B.imag);						      \
    }

#define	C_ADDEQ(A,B,C) {						      \
	(A).real = (B.real) + (C.real);					      \
	(A).imag = (B.imag) + (C.imag);					      \
    }

#define	C_SUB(A,B) {							      \
	(A).real -= (B.real);						      \
	(A).imag -= (B.imag);						      \
    }

#define	C_SUBEQ(A,B,C) {						      \
	(A).real = (B.real) - (C.real);					      \
	(A).imag = (B.imag) - (C.imag);					      \
    }



#endif /*CMPLX*/
