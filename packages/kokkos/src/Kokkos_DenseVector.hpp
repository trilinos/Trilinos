//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DENSEVECTOR_H
#define KOKKOS_DENSEVECTOR_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_Vector.hpp"


namespace Kokkos {

//! Kokkos::DenseVector: Kokkos vector base class.

/*! The Kokkos::DenseVector specifies the interface that any vector class interfacing to the Kokkos 
  Operators classes must implement.

  At this time, the primary function provided by Kokkos::DenseVector is access to vector data.

*/    

  template<typename OrdinalType, typename ScalarType>
  class DenseVector: public virtual Vector<OrdinalType,ScalarType> {
  public:

    //@{ \name Constructors/Destructor.

    //! Default constructor
    DenseVector(void):
      dataInitialized_(false),
      length_(0),
      inc_(0),
      values_(0) {};
  
    //! Copy constructor.
    DenseVector(const DenseVector& source):
      dataInitialized_(source.dataInitialized_),
      length_(source.length_),
      inc_(source.inc_),
      values_(source.values_) {};

    //! DenseVector Destructor
    virtual ~DenseVector(){};
    //@}

    //@{ \name Initialization methods
	
    //! Initialize using a pointer
    /*!
      This is the only way to initialize a Kokkos::DenseVector object.
      \param length (In)  Length of vector.
      \param values (In)  Pointer to values.
      \param inc (In) The increment between two elements in the vector.  
                         Typically this value should be set to 1, which is the default value.

      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(OrdinalType length, ScalarType * values, OrdinalType inc = 1) {
      length_ = length;
      inc_ = inc;
      values_ = values;
      dataInitialized_ = true;
      return(0);
      };
	
    //@}

    //@{ \name DenseVector access methods.

    //! Returns a pointer to an array of values in the vector.
    /*! Extract a pointer to the values in the vector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the vector object itself.  The getInc() method 
	should be used to access values, especially if getInc() != 1.
    */
    virtual ScalarType * getValues() const {return(values_);};
	
	
    //@}

    //@{ \name DenseVector Attribute access methods.
	
    //! Length of vector
    virtual OrdinalType getLength() const {return(length_);};
	
    //! Increment between entries in the vector, normally = 1.
    virtual OrdinalType getInc() const {return(inc_);};
	
    //@}

    bool dataInitialized_;
    OrdinalType length_;
    OrdinalType inc_;

    ScalarType * values_;
  };

} // namespace Kokkos
#endif /* KOKKOS_DENSEVECTOR_H */
