//@HEADER
// ************************************************************************
// 
//                 Kokkos: A Fast Kernel Package
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

#ifndef KOKKOS_OSKIVECTOR_H
#define KOKKOS_OSKIVECTOR_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_OskiMultiVector.hpp"


namespace Kokkos {

//! Kokkos::OskiVector: Implementation of the abstract Kokkos:Vector that wraps an oski_vecview_t object with one vector.  The class is derived from the Kokkos:OskiMultiVector class.  As such it has full access to all of the functionality provided by the Kokkos:OskiMultiVector class. 

// Ask Mike about whether or not it would be better to have this implement
// the OskiVector base class.

/*! 
  At this time, the primary function provided by Kokkos::DenseVector is wrapping an oski_vecview_t object and providing access to its entries.

*/    

  template<typename OrdinalType, typename ScalarType>
//  class OskiVector: public virtual Vector<OrdinalType,ScalarType> {
    class OskiVector: public OskiMultiVector<OrdinalType,ScalarType> {
  public:

    //! @name Constructors/Destructor

    //@{

    //! Default constructor
    OskiVector(void):
      OskiMultiVector<OrdinalType,ScalarType>() {};
  
    //! Copy constructor.
    OskiVector(const OskiVector& source):
      OskiMultiVector<OrdinalType, ScalarType>(source) {};

    // Based on the Epetra_Vector class, it doesn't look like I need to call
    // the OskiMultiVector Destructor
    //! OskiVector Destructor
    virtual ~OskiVector(){};
    //@}

    //! @name Initialization methods

    //@{

    //! Initialize using a pointer
    /*!
      This is the only way to initialize a Kokkos::OskiVector object.
      \param length (In)  Length of vector.
      \param values (In)  Pointer to values.
      \param inc (In) The increment between two elements in the vector.  
                         Typically this value should be set to 1, which is the default value.

      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(OrdinalType length, ScalarType * values, OrdinalType inc = 1) {
     OskiMultiVector<OrdinalType,ScalarType>::initializeValues(length,1,values,inc); 
     return(0);
      };
	
    //@}

    //! @name OskiVector access methods

    //@{

    //! Returns a pointer to an array of values in the vector.
    /*! Extract a pointer to the values in the vector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the vector object itself.  The getInc() method 
	should be used to access values, especially if getInc() != 1.
    */
//consider changing this to *.  Oski uses *, the ** is just to stay consistant
// with MultiVector.  That would have to be changed as well.
    virtual ScalarType ** getValues() const {return(OskiMultiVector<OrdinalType,ScalarType>::getValues());};
	
	
    //@}

    //! @name OskiVector Attribute access methods

    //@{
	
    //! Length of vector
    virtual OrdinalType getLength() const {return(OskiMultiVector<OrdinalType,ScalarType>::getNumRows());};
	
    //! Increment between entries in the vector, normally = 1.
    virtual OrdinalType getInc() const {return(OskiMultiVector<OrdinalType,ScalarType>::getRowInc());};

    //! Underlying OSKI Vector
        virtual oski_vecview_t getX_view() const{return(OskiMultiVector<OrdinalType,ScalarType>::getX_view());};

    //@}

//  protected:
//    oski_vecview_t x_view_;

//    bool dataInitialized_;
//    OrdinalType length_;
//    OrdinalType inc_;

//    ScalarType * values_;
  };

} // namespace Kokkos
#endif /* KOKKOS_OSKIVECTOR_H */

