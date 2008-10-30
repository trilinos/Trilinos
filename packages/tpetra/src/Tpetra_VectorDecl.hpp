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

#include <Teuchos_Object.hpp>
#include <Teuchos_ArrayView.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of VectorData, needed to prevent circular inclusions
  template<typename Ordinal, typename Scalar> class VectorData;
#endif

  //! Tpetra::Vector: A class for constructing and using sparse vectors.

  /*! Vector is templated on Scalar for the vector entries, and on Ordinal 
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

  template<typename Ordinal, typename Scalar>
  class Vector : public DistObject<Ordinal,Scalar> {

  public:
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    Vector(const Map<Ordinal> &map);
  
    //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
    Vector(const Teuchos::ArrayView<const Scalar> &values, const Map<Ordinal> & map);

    //! Copy constructor.
    Vector(const Vector<Ordinal,Scalar> &source);

    //! Destructor.  
    ~Vector();

    //@}

    //! @name Post-Construction Modification Routines
    //@{ 

    //! Submit entries. Values submitted will be summed with existing values.
    void submitEntries(const Teuchos::ArrayView<const Ordinal> &indices,
                       const Teuchos::ArrayView<const Scalar>  &values);

    //! Set all entries to scalarValue.
    void setAllToScalar(const Scalar &value);

    //! Set all entries to random values.
    void setAllToRandom();

    //@}


    //! @name Extraction Methods
    //@{ 

    //! Put vector entries into user array (copy)
    // void extractCopy(Scalar* userArray) const;

    //! Put pointers to vector entries into user array (view)
    // void extractView(Scalar** userPointerArray) const;

    //@}


    //! @name Mathematical Methods
    //@{ 

    //! Returns result of dot product, \e result = this.x
    Scalar dotProduct(const Vector<Ordinal, Scalar> &x) const;

    //! Changes this vector to elementwise absolute values of x.
    void absoluteValue(const Vector<Ordinal,Scalar> &x);

    //! Changes this vector to element-wise reciprocal values of x.
    void reciprocal(const Vector<Ordinal,Scalar> &x);

    //! Scale the current values of a vector, \e this = scalarThis*\e this.
    void scale(const Scalar &scalarThis);

    //! Replace vector values with scaled values of x, \e this = scalarX*x.
    void scale(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x);

    //! Update vector values with scaled values of x, \e this = scalarThis*\e this + scalarX*x.
    void update(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x, const Scalar &scalarThis);

    //! Update vector with scaled values of x and y, \e this = scalarThis*\e this + scalarX*x + scalarY*y.
    void update(const Scalar &scalarX, const Vector<Ordinal,Scalar> &x, 
                const Scalar &scalarY, const Vector<Ordinal,Scalar> &y, 
                const Scalar &scalarThis);

    //! Compute 1-norm of vector.
    Scalar norm1() const;

    //! Compute 2-norm of vector.
    Scalar norm2() const;

    //! Compute Infinity-norm of vector.
    Scalar normInf() const;

    //! Compute Weighted 2-norm (RMS Norm) of vector.
    Scalar normWeighted(const Vector<Ordinal,Scalar> &weights) const;

    //! Compute minimum value of vector.
    Scalar minValue() const;

    //! Compute maximum value of vector.
    Scalar maxValue() const;

    //! Compute mean (average) value of vector.
    Scalar meanValue() const;

    //! Vector multiplication (elementwise) 
    /*! \e this = scalarThis*\e this + scalarXY*x@y, where @ represents elementwise multiplication. */
    void elementwiseMultiply(const Scalar &scalarXY, const Vector<Ordinal,Scalar> &x, const Vector<Ordinal,Scalar> &y, 
                             const Scalar &scalarThis);

    //! Reciprocal multiply (elementwise)
    /*! \e this = scalarThis*\e this + scalarXY*y@x, where @ represents elementwise division. */
    void elementwiseReciprocalMultiply(Scalar scalarXY, const Vector<Ordinal, Scalar> &x, const Vector<Ordinal, Scalar> &y, 
                                       const Scalar &scalarThis);

    //@}


    //! @name Random number utilities
    //@{ 

    //! Get seed
    const Scalar & getSeed() const;

    //! Set seed
    void setSeed(const Scalar &seed);

    //@}


    //! @name Element access methods
    //@{ 

    //! [] operator, nonconst version
    Scalar& operator[](Ordinal index);

    //! [] operator, const version
    const Scalar & operator[](Ordinal index) const;

    //@}


    //! @name Attribute access methods
    //@{ 

    //! Returns number of vector entries owned by this image.
    Ordinal getNumMyEntries() const;

    //! Returns number of vector entries across all images.
    Ordinal getNumGlobalEntries() const;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method, used by overloaded << operator.
    void print(std::ostream &os) const;

    void printValues(std::ostream &os) const;

    //@}

    //! @name Misc. 
    //@{ 

    //! Returns a const reference to the VectorSpace this Vector belongs to.
    const Map<Ordinal> & getMap() const;

    //! Assignment Operator
    Vector<Ordinal,Scalar> & operator=(const Vector<Ordinal,Scalar> &source);

    //@}

    //! @name Expert/Developer Use Only.
    //@{ 

    // Returns pointer to Scalar array inside of scalarArray
    Teuchos::ArrayView<Scalar> scalarPointer();

    Teuchos::ArrayView<const Scalar> scalarPointer() const;

    //@}

  private:

    Teuchos::RCP<VectorData<Ordinal,Scalar> > VectorData_;

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Ordinal,Scalar> & sourceObj);

    void copyAndPermute(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numSameIDs,
               Ordinal numPermuteIDs,
               const Teuchos::ArrayView<const Ordinal> & permuteToLIDs,
               const Teuchos::ArrayView<const Ordinal> & permuteFromLIDs);

    void packAndPrepare(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numExportIDs,
               const Teuchos::ArrayView<const Ordinal> & exportLIDs,
               const Teuchos::ArrayView<Scalar> & exports,
               Ordinal &packetSize,
               Distributor<Ordinal> &distor);
  
    void unpackAndCombine(Ordinal numImportIDs,
               const Teuchos::ArrayView<const Ordinal> & importLIDs,
               const Teuchos::ArrayView<const Scalar> & imports,
               Distributor<Ordinal> &distor,
               CombineMode CM);

  }; // class Vector

} // namespace Tpetra

#endif // TPETRA_VECTOR_DECL_HPP
