/*Paul
19-Nov-2002 PacketTraits written, based on ScalarTraits form
*/

#ifndef _TPETRA_PACKETTRAITS_HPP_
#define _TPETRA_PACKETTRAITS_HPP_

namespace Tpetra {
  /** The Tpetra PacketTraits file.
			
	For the most general type 'class T', we define aborting functions, 
	which should restrict implementations from using traits other than the
	specializations defined below.
  */
	
	template<class T>
	struct PacketTraits {
		static inline int unsupportedType() {
#ifndef TPETRA_NO_ERROR_REPORTS
			cerr << endl << "Tpetra::PacketTraits: unsupported packet type." << endl;
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

#ifndef NO_COMPLEX
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

#endif // NO_COMPLEX

} // namespace Tpetra

#endif // _TPETRA_ORDINALTRAITS_HPP_
