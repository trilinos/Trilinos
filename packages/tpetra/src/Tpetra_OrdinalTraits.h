/*Paul
19-Nov-2002 OrdinalTraits written, based on ScalarTraits form
*/

#ifndef _TPETRA_ORDINALTRAITS_H_
#define _TPETRA_ORDINALTRAITS_H_

namespace Tpetra {
  /** The Tpetra OrdinalTraits file.
			
	For the most general type 'class T', we define aborting functions, 
	which should restrict implementations from using traits other than the
	specializations defined below.
  */
	
	template<class T>
	struct OrdinalTraits {

		static inline int unsupportedType() {
#ifndef TPETRA_NO_ERROR_REPORTS
			cerr << endl << "Tpetra::OrdinalTraits: unsupported ordinal type." << endl;
#endif
			return(-1);
		}

		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if ordinal traits machine parameters defined
		static inline T zero()                     {throw(unsupportedType());};
		static inline T one()                      {throw(unsupportedType());};
		static inline const char* name()           {throw(unsupportedType());};
	};
	
	template<>
	struct OrdinalTraits<int> {

		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if scalar traits machine parameters defined 
		static inline int zero()                   {return(0);};
		static inline int one()                    {return(1);};
		static inline const char* name()           {return("int");};
	};

} // namespace Tpetra

#endif // _TPETRA_ORDINALTRAITS_H_
