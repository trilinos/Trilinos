#ifndef _TEUCHOS_PACKETTRAITS_HPP_
#define _TEUCHOS_PACKETTRAITS_HPP_

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

namespace Teuchos {
  /** The Teuchos PacketTraits file.
			
	For the most general type 'class T', we define aborting functions, 
	which should restrict implementations from using traits other than the
	specializations defined below.
  */
	
	template<class T>
	struct PacketTraits {
		static inline int unsupportedType() {
#ifndef TEUCHOS_NO_ERROR_REPORTS
			cerr << endl << "Teuchos::PacketTraits: unsupported packet type." << endl;
#endif
			return(-1);
		}
		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if packet traits machine parameters defined
		static inline int packetSize()             {throw(unsupportedType());};
		static inline const char* name()           {throw(unsupportedType());};
	};
	
	template<>
	struct OrdinalTraits<int> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(int));};
		static inline const char* name()           {return("int");};
	};

	template<>
	struct OrdinalTraits<float> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(float));};
		static inline const char* name()           {return("float");};
	};

	template<>
	struct OrdinalTraits<double> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(double));};
		static inline const char* name()           {return("double");};
	};

#if (defined(HAVE_COMPLEX) || defined(HAVE_COMPLEX_H)) && defined(HAVE_TEUCHOS_EXPERIMENTAL)
#include <complex>

	template<>
	struct OrdinalTraits<complex<float> > {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(complex<float>));};
		static inline const char* name()           {return("complex<float>");};
	};

	template<>
	struct OrdinalTraits<complex<double> > {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(complex<double>));};
		static inline const char* name()           {return("complex<double>");};
	};

#endif // HAVE_COMPLEX || HAVE_COMPLEX_H

} // namespace Teuchos

#endif // _TEUCHOS_PACKETTRAITS_HPP_
