#ifndef _ESI_scalarTraits_h_
#define _ESI_scalarTraits_h_

#include <math.h> //for the fabs function...

namespace esi {
/** The ESI scalarTraits file.

    For the most general type 'class T', we define aborting functions, 
		which should restrict implementations from using traits other than the
		specializations defined below.
*/

template<class T>
struct scalarTraits {

	typedef T magnitude_type;

	static inline magnitude_type magnitude(T a) { 
		cout << "esi::scalarTraits: unsupported scalar type." << endl; 
		abort(); 
		return(a);
	};

	static inline T random() {
		cout << "esi::scalarTraits: unsupported scalar type." << endl; 
		abort(); 
		return(a);
	};
	
	static inline const char * name() {
		cout << "esi::scalarTraits: unsupported scalar type." << endl; 
		abort(); 
		return(0); 
	};
};

template<>
struct scalarTraits<esi_real4> {

	typedef esi_real4 magnitude_type;
	
	static inline magnitude_type magnitude(esi_real4 a) { 
		return(fabs(a)); 
	};

	static inline esi_real4 random() {
		esi_real4 rnd = (esi_real4)rand()/RAND_MAX;
		return( (esi_real4)(-1.0 + 2.0*rnd) );
	};

	static inline const char * name() { 
		return("esi_real4"); 
	};
};


template<>
struct scalarTraits<esi_real8> {

	typedef esi_real8 magnitude_type;

	static magnitude_type magnitude(esi_real8 a) {
		return(fabs(a));
	};

	static inline esi_real8 random() {
		esi_real8 rnd = (esi_real8)rand()/RAND_MAX;
		return( (esi_real8)(-1.0 + 2.0*rnd) );
	};

	static inline const char * name() {
		return("esi_real8");
	};
};

/*
  If we're using a Sun compiler, version earlier than 5.0,
  then complex isn't available.
*/
#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define NO_COMPLEX
#endif

/*
  If we're using the tflops Portland Group compiler, then complex isn't
  available. (As of July 21, 2000. abw)
*/
#if defined(__PGI) && defined(__i386)
#define NO_COMPLEX
#endif

#define NO_COMPLEX
#ifndef NO_COMPLEX

#include <complex>

template<>
struct scalarTraits< std::complex<esi_real4> > {

	typedef esi_real4 magnitude_type;

	static magnitude_type magnitude(std::complex<esi_real4> a) {
		return(std::abs(a));
	};

	static inline std::complex<esi_real4> random() {
		esi_real4 rnd1 = (esi_real4)rand()/RAND_MAX;
		esi_real4 rnd2 = (esi_real4)rand()/RAND_MAX;
		return( std::complex<esi_real4>(-1.0+2.0*rnd1, -1.0+2.0*rnd2) );
	};

	static inline const char * name() {
		return("complex<esi_real4>");
	};
};

template<>
struct scalarTraits< std::complex<esi_real8> > {

	typedef esi_real8 magnitude_type;

	static magnitude_type magnitude(std::complex<esi_real8> a) {
		return(std::abs(a));
	};
	
	static inline std::complex<esi_real8> random() {
		esi_real8 rnd1 = (esi_real8)rand()/RAND_MAX;
		esi_real8 rnd2 = (esi_real8)rand()/RAND_MAX;
		return( std::complex<esi_real8>(-1.0+2.0*rnd1, -1.0+2.0*rnd2) );
	};

	static inline const char * name() {
		return("complex<esi_real8>");
	};
};

#endif  /* NO_COMPLEX */


/*
  If we're using a Sun compiler, version earlier than 5.0,
  then 'typename' isn't available.
*/
#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define TYPENAME
#else
#define TYPENAME typename
#endif

};     // esi namespace
#endif

