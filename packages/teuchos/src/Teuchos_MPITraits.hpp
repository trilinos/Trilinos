#ifndef TEUCHOS_MPITRAITS_H
#define TEUCHOS_MPITRAITS_H

#include "Teuchos_MPIComm.hpp"

namespace Teuchos
{
	using std::string;

	/** \ingroup MPI 
	 * MPITraits let us bind MPI types to C++ types.
	 */
	template <class T> class MPITraits
		{
		public:
			/** return the MPI data type of the template argument */
			static int type();
		};
	
	/** \ingroup MPI 
	 * Binds MPI_INT to int
	 */
	template <> class MPITraits<int>
		{
		public:
			/** return the MPI data type of the template argument */
			static int type() {return MPIComm::INT;}
		};
	
	/** \ingroup MPI 
	 * Binds MPI_FLOAT to float
	 */
	template <> class MPITraits<float>
		{
		public:
			/** return the MPI data type of the template argument */
			static int type() {return MPIComm::FLOAT;}
		};
	
	/** \ingroup MPI 
	 * Binds MPI_DOUBLE to double
	 */
	template <> class MPITraits<double>
		{
		public:
			/** return the MPI data type of the template argument */
			static int type() {return MPIComm::DOUBLE;}
		};
	
	/** \ingroup MPI 
	 * Binds MPI_CHAR to char
	 */
	template <> class MPITraits<char>
		{
		public:
			/** return the MPI data type of the template argument */
			static int type() {return MPIComm::CHAR;}
		};
	
}

#endif
