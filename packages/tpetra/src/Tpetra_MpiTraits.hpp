// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_MPITRAITS_HPP
#define TPETRA_MPITRAITS_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <mpi.h>

namespace Tpetra {
  /** The Tpetra MPITraits file.
			
	This is a traits file used by MpiComm. It provides traits
	specializations for the built-in types that MPI can handle.
  */	

  template<class T>
  struct UndefinedScalarTraits {
    //! This function should not compile if there is an attempt to instantiate!
    static inline T notDefined() { return T::this_type_is_missing_a_specialization(); };
  };
  
	template<class T>
	struct MpiTraits {
    	static inline MPI_Datatype dataType()    {return UndefinedScalarTraits<T>::notDefined();};
	};
	
	template<>
	struct MpiTraits<int> {
    	static inline MPI_Datatype dataType()    {return(MPI_INT);}; 
	};

	template<>
	struct MpiTraits<float> {
    	static inline MPI_Datatype dataType()    {return(MPI_FLOAT);}; 
	};

	template<>
	struct MpiTraits<double> {
    	static inline MPI_Datatype dataType()    {return(MPI_DOUBLE);}; 
	};
   
} // namespace Tpetra

#endif // TPETRA_MPITRAITS_HPP
