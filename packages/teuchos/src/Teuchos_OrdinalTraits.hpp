/*Paul
19-Nov-2002 OrdinalTraits written, based on ScalarTraits form
*/

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_ORDINALTRAITS_HPP_
#define _TEUCHOS_ORDINALTRAITS_HPP_

namespace Teuchos {
  /** The Teuchos OrdinalTraits file.
			
	For the most general type 'class T', we define aborting functions, 
	which should restrict implementations from using traits other than the
	specializations defined below.
  */
	
	template<class T>
	struct OrdinalTraits {

		static inline int unsupportedType() {
#ifndef TEUCHOS_NO_ERROR_REPORTS
			cerr << endl << "Teuchos::OrdinalTraits: unsupported ordinal type." << endl;
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

		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if ordinal traits machine parameters defined 
		static inline int zero()                   {return(0);};
		static inline int one()                    {return(1);};
		static inline const char* name()           {return("int");};
	};

} // namespace Teuchos

#endif // _TEUCHOS_ORDINALTRAITS_HPP_
