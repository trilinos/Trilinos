// ////////////////////////////////////////////////////////////////////////
// Teuchos_RawMPITraits.hpp

#ifndef TEUCHOS_RAW_MPI_TRAITS_H
#define TEUCHOS_RAW_MPI_TRAITS_H

#include "Teuchos_ConfigDefs.hpp"

/** \file Teuchos_RawMPITraits.hpp
 *
 * Note, this file should only be included if MPI is available is defined.
 */

namespace Teuchos {

///
/** Declaration of a traits class that returns raw data types.
 *
 * Note that this class should not compile if it is instantiated by
 * accident.
 */
template <class T> class RawMPITraits {
public:
	///
	static MPI_Datatype type() { bool *junk1; T *junk2 = &junk1; return MPI_DATATYPE_NULL; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>char</tt>
 */
template <> class RawMPITraits<char> {
public:
	///
	static MPI_Datatype type() { return MPI_CHAR; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>int</tt>
 */
template <> class RawMPITraits<int> {
public:
	///
	static MPI_Datatype type() { return MPI_INT; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>float</tt>
 */
template <> class RawMPITraits<float> {
public:
	///
	static MPI_Datatype type() { return MPI_FLOAT; }
};

///
/** Specialization of <tt>RawMPITraits</tt> for <tt>double</tt>
 */
template <> class RawMPITraits<double> {
public:
	///
	static MPI_Datatype type() { return MPI_DOUBLE; }
};

} // namespace Teuchos

#endif // TEUCHOS_RAW_MPI_TRAITS_H
