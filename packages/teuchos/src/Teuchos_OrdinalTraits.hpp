// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_ORDINALTRAITS_HPP_
#define _TEUCHOS_ORDINALTRAITS_HPP_

/*! \file Teuchos_OrdinalTraits.hpp
   \brief Defines basic traits for the ordinal field type 
*/

/*! \struct Teuchos::OrdinalTraits
    \brief This structure defines some basic traits for the ordinal field type.

    Ordinal traits are an essential part of templated codes.  This structure offers
    the basic traits of the templated ordinal type, like defining zero and one.

    For the general type, or default implementation, an aborting function
    is defined which should restrict implementations from using ordinal traits other than
    the defined specializations.

    \note The defined specializations for OrdinalTraits are: \c int and \c long \c int.
*/

namespace Teuchos {
	
	template<class T>
	struct OrdinalTraits {

		//! Aborting function to restrict non-supported implementations of OrdinalTraits.
		static inline int unsupportedType() {
#ifndef TEUCHOS_NO_ERROR_REPORTS
			cerr << endl << "Teuchos::OrdinalTraits: unsupported ordinal type." << endl;
#endif
			return(-1);
		}

		//! Allows testing to see if ordinal traits machine parameters are defined.
		static inline bool haveMachineParameters() {return(false);}; 
		
		//! Returns representation of zero for this ordinal type.
		static inline T zero()                     {throw(unsupportedType());};

		//! Returns representation of one for this ordinal type.
		static inline T one()                      {throw(unsupportedType());};

		//! Returns name of this ordinal type.
		static inline const char* name()           {throw(unsupportedType());};
	};
	
#ifndef DOXYGEN_SHOULD_SKIP_THIS

	template<>
	struct OrdinalTraits<int> {

		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if ordinal traits machine parameters defined 
		static inline int zero()                   {return(0);};
		static inline int one()                    {return(1);};
		static inline const char* name()           {return("int");};
	};

	template<>
	struct OrdinalTraits<long int> {

		static inline bool haveMachineParameters() {return(false);}; // Allows testing to see if ordinal traits machine parameters defined 
		static inline long int zero()                   {return((long int)0);};
		static inline long int one()                    {return((long int)1);};
		static inline const char* name()           {return("long int");};
	};

#endif

} // namespace Teuchos

#endif // _TEUCHOS_ORDINALTRAITS_HPP_
