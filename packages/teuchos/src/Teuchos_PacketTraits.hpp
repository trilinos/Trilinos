#ifndef _TPETRA_PACKETTRAITS_HPP_
#define _TPETRA_PACKETTRAITS_HPP_

// Kris
// 07.08.03 -- Move into Teuchos package/namespace
// 01.16.04 -- Move back into Tpetra package/namespace per discussion with Mike Heroux.

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {
  /** The Tpetra PacketTraits file.
			
	For the most general type 'class T', we define aborting functions, 
	which should restrict implementations from using traits other than the
	specializations defined below.
  */	

	template<class T>
	struct PacketTraits {
		static inline int unsupportedType() {
#ifndef TEUCHOS_NO_ERROR_REPORTS
			cerr << endl << "Tpetra::PacketTraits: unsupported packet type." << endl;
#endif
			return(-1);
		}
		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if packet traits machine parameters defined
		static inline int packetSize()             {throw(unsupportedType());};
		static inline const char* name()           {throw(unsupportedType());};
	};
	
	template<>
	struct PacketTraits<int> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(int));};
		static inline const char* name()           {return("int");};
	};

	template<>
	struct PacketTraits<float> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(float));};
		static inline const char* name()           {return("float");};
	};

	template<>
	struct PacketTraits<double> {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(double));};
		static inline const char* name()           {return("double");};
	};

#if (defined(HAVE_COMPLEX) || defined(HAVE_COMPLEX_H))

	template<>
	struct PacketTraits<complex<float> > {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(complex<float>));};
		static inline const char* name()           {return("complex<float>");};
	};

	template<>
	struct PacketTraits<complex<double> > {
		static inline bool haveMachineParameters() {return(true);}; // Allows testing to see if packet traits machine parameters defined 
		static inline int packetSize()             {return(sizeof(complex<double>));};
		static inline const char* name()           {return("complex<double>");};
	};

#endif // HAVE_COMPLEX || HAVE_COMPLEX_H

} // namespace Tpetra

#endif // _TPETRA_PACKETTRAITS_HPP_
