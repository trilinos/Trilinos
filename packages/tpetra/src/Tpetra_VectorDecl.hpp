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

#ifndef TPETRA_VECTOR_DECL_HPP
#define TPETRA_VECTOR_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_OutputObject.hpp" // for STL vector, algorithm, numeric
#include "Tpetra_VectorSpace.hpp"
#include <Teuchos_CompObject.hpp>
#include <Teuchos_BLAS.hpp>

namespace Tpetra {

  // forward declaration of VectorData, needed to prevent circular inclusions
  template<typename OrdinalType, typename ScalarType> class VectorData;

  //! Tpetra::Vector: A class for constructing and using sparse vectors.

  /*! Vector is templated on ScalarType for the vector entries, and on OrdinalType 
      for the vector indices. A VectorSpace object is needed for all Vector objects.

    Vector entries can only be accessed through their local index values. 
    Global index values can be converted to local indices by using the 
    VectorSpace::getLocalIndex method.

    Note that for most of the mathematical methods that set \e this to the result of
    an operation on vectors passed as parameters, the \e this vector can be used 
    as one of the parameters (unless otherwise specified).
    
    Vector error codes (positive for non-fatal, negative for fatal):
    <ol>
    <li> +1  Specified vector index not found on this image.
    <li> +2  Vector sizes do not match.
    <li> +3  Cannot perform that operation on an empty vector.
    <li> -1  Invalid number of entries passed to constructor.
    <li> -99 Internal Vector error. Contact developer.
    </ol>
  */
  
  template<typename OrdinalType, typename ScalarType>
  class Vector : public Teuchos::CompObject, public DistObject<OrdinalType, ScalarType>, public virtual OutputObject {

  public:
  
    //@{ \name Constructor/Destructor Methods

    //! Sets all vector entries to zero.
    Vector(VectorSpace<OrdinalType, ScalarType> const& VectorSpace);
  
    //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
    Vector(const ScalarType* vectorEntries, OrdinalType numEntries, VectorSpace<OrdinalType, ScalarType> const& VectorSpace);

    //! Copy constructor.
    Vector(Vector<OrdinalType, ScalarType> const& Source);

    //! Destructor.  
    ~Vector();

    //@}

    //@{ \name Post-Construction Modification Routines

    //! Submit entries. Values submitted will be summed with existing values.
    void submitEntries(OrdinalType numEntries, OrdinalType const* indices, ScalarType const* values);

    //! Set all entries to scalarValue.
    void setAllToScalar(ScalarType const value);

    //! Set all entries to random values.
    void setAllToRandom();

    //@}


    //@{ \name Extraction Methods

    //! Put vector entries into user array (copy)
    void extractCopy(ScalarType* userArray) const;

    //! Put pointers to vector entries into user array (view)
    void extractView(ScalarType** userPointerArray) const;

    //@}


    //@{ \name Mathematical Methods

    //! Returns result of dot product, \e result = this.x
    ScalarType dotProduct(Vector<OrdinalType, ScalarType> const& x) const;

    //! Changes this vector to elementwise absolute values of x.
    void absoluteValue(Vector<OrdinalType, ScalarType> const& x);

    //! Changes this vector to element-wise reciprocal values of x.
    void reciprocal(Vector<OrdinalType, ScalarType> const& x);

    //! Scale the current values of a vector, \e this = scalarThis*\e this.
    void scale(ScalarType scalarThis);

    //! Replace vector values with scaled values of x, \e this = scalarX*x.
    void scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x);

    //! Update vector values with scaled values of x, \e this = scalarThis*\e this + scalarX*x.
    void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarThis);

    //! Update vector with scaled values of x and y, \e this = scalarThis*\e this + scalarX*x + scalarY*y.
    void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarY, 
          Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis);

    //! Compute 1-norm of vector.
    ScalarType norm1() const;

    //! Compute 2-norm of vector.
    ScalarType norm2() const;

    //! Compute Infinity-norm of vector.
    ScalarType normInf() const;

    //! Compute Weighted 2-norm (RMS Norm) of vector.
    ScalarType normWeighted(Vector<OrdinalType, ScalarType> const& weights) const;

    //! Compute minimum value of vector.
    ScalarType minValue() const;

    //! Compute maximum value of vector.
    ScalarType maxValue() const;

    //! Compute mean (average) value of vector.
    ScalarType meanValue() const;

    //! Vector multiplication (elementwise) 
    /*! \e this = scalarThis*\e this + scalarXY*x@y, where @ represents elementwise multiplication. */
    void elementwiseMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
                 Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis);

    //! Reciprocal multiply (elementwise)
    /*! \e this = scalarThis*\e this + scalarXY*y@x, where @ represents elementwise division. */
    void elementwiseReciprocalMultiply(ScalarType scalarXY, 
                       Vector<OrdinalType, ScalarType> const& x, 
                       Vector<OrdinalType, ScalarType> const& y, 
                       ScalarType scalarThis);

    //@}


    //@{ \name Random number utilities

    //! Get seed
    ScalarType getSeed() const;

    //! Set seed
    void setSeed(ScalarType seed);

    //@}


    //@{ \name Element access methods

    //! [] operator, nonconst version
    ScalarType& operator[](OrdinalType index);

    //! [] operator, const version
    ScalarType const& operator[](OrdinalType index) const;

    //@}


    //@{ \name Attribute access methods

    //! Returns number of vector entries owned by this image.
    OrdinalType getNumMyEntries() const;

    //! Returns number of vector entries across all images.
    OrdinalType getNumGlobalEntries() const;

    //@}


    //@{ \name I/O methods

    //! Print method, used by overloaded << operator.
    void print(ostream& os) const;
  
    void printValues(ostream& os) const;

    //@}

    //@{ \name Misc. 

    //! Returns a const reference to the VectorSpace this Vector belongs to.
    VectorSpace<OrdinalType, ScalarType> const& vectorSpace() const;

    //! Assignment Operator
    Vector<OrdinalType, ScalarType>& operator = (Vector<OrdinalType, ScalarType> const& Source);

    //@}

    //@{ \name Expert/Developer Use Only.

    // Returns pointer to ScalarType array inside of scalarArray
    ScalarType* scalarPointer();

    ScalarType const* scalarPointer() const;

    //@}

  private:

    // Accessor for BLAS
    Teuchos::BLAS<OrdinalType, ScalarType> const& BLAS() const;

    // Accessors for scalarArray
    std::vector<ScalarType>& scalarArray();
    std::vector<ScalarType>const & scalarArray() const;

    Teuchos::RCP< VectorData<OrdinalType, ScalarType> > VectorData_;

    // four functions needed for DistObject derivation
    bool checkSizes(DistObject<OrdinalType, ScalarType> const& sourceObj);

    int copyAndPermute(DistObject<OrdinalType, ScalarType> const& sourceObj,
               OrdinalType const numSameIDs,
               OrdinalType const numPermuteIDs,
               std::vector<OrdinalType> const& permuteToLIDs,
               std::vector<OrdinalType> const& permuteFromLIDs);

    int packAndPrepare(DistObject<OrdinalType, ScalarType> const& sourceObj,
               OrdinalType const numExportIDs,
               std::vector<OrdinalType> const& exportLIDs,
               std::vector<ScalarType>& exports,
               OrdinalType& packetSize,
               Distributor<OrdinalType> const& distor);
  
    int unpackAndCombine(OrdinalType const numImportIDs,
               std::vector<OrdinalType> const& importLIDs,
               std::vector<ScalarType> const& imports,
               Distributor<OrdinalType> const& distor,
               CombineMode const CM);

  }; // class Vector

} // namespace Tpetra

#endif // TPETRA_VECTOR_DECL_HPP
