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

#ifndef TPETRA_PACKETTRAITS_HPP
#define TPETRA_PACKETTRAITS_HPP

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {
  /*! The Tpetra PacketTraits file.
			
    For the general type, or default implementation, an aborting function
    is defined which should restrict implementations from using ordinal traits other than
    the defined specializations.
  */	

  template<typename T>
  struct UndefinedPacketTraits {
    //! This function should not compile if there is an attempt to instantiate!
    static inline T notDefined() { return T::this_type_is_missing_a_specialization(); };
  };
  
	template<typename T>
	struct PacketTraits {
		static inline int packetSize()              {return UndefinedPacketTraits<T>::notDefined();};
		static inline std::string name()            {return UndefinedPacketTraits<T>::notDefined();};
	};
	
#ifndef DOXYGEN_SHOULD_SKIP_THIS

	template<>
	struct PacketTraits<int> {
		static inline int packetSize()              {return(sizeof(int));};
		static inline std::string name()            {return("int");};
	};

	template<>
	struct PacketTraits<unsigned int> {
		static inline int packetSize()              {return(sizeof(unsigned int));};
		static inline std::string name()            {return("unsigned int");};
	};

	template<>
	struct PacketTraits<long int> {
		static inline int packetSize()              {return(sizeof(long int));};
		static inline std::string name()            {return("long int");};
	};

	template<>
	struct PacketTraits<unsigned long int> {
		static inline int packetSize()              {return(sizeof(unsigned long int));};
		static inline std::string name()            {return("unsigned long int");};
	};

	template<>
	struct PacketTraits<float> {
		static inline int packetSize()              {return(sizeof(float));};
		static inline std::string name()            {return("float");};
	};

	template<>
	struct PacketTraits<double> {
		static inline int packetSize()              {return(sizeof(double));};
		static inline std::string name()            {return("double");};
	};

	template<>
	struct PacketTraits<long double> {
		static inline int packetSize()              {return(sizeof(long double));};
		static inline std::string name()            {return("long double");};
	};

  template<typename T>
	struct PacketTraits< complex<T> > {
		static inline int packetSize()              {return(sizeof(complex<T>));};
		static inline std::string name()            {return("complex<" + PacketTraits<T>::name() + ">");};
	};
  

#endif // DOXYGEN_SHOULD_SKIP_THIS
   
} // namespace Tpetra

#endif // TPETRA_PACKETTRAITS_HPP
