// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_VECTOR_H
#define NOX_EPETRA_VECTOR_H

#include "NOX_Abstract_Vector.H"     // base class
#include "Teuchos_RCP.hpp"   // class data element
#include "NOX_Epetra_VectorSpace.H"  // class data element

// Forward declarations
class Epetra_Vector;
namespace NOX {
  namespace Epetra {
    class VectorSpace;
  }
}

namespace NOX {

//! %NOX %Epetra support.
namespace Epetra {

//! Implementation of NOX::Abstract::Vector for %Epetra vectors.
class Vector : public virtual NOX::Abstract::Vector {

 public:

  //! Type of memory management to use when constructing the vector.
  enum MemoryType {
    //! Keeps a pointer to and uses the actual Epetra_Vector passed in.
    CreateView,
    //! Allocates a new underlying Epetra_Vector object.
    CreateCopy
  };

  /*! \brief Constructor that creates a COPY or VIEW of the Epetra_Vector.

    NOTE: This ctor should just always create a view.  It should be
    implicit from the fact that a RCP object is being passed
    in that a persisting relationship is present.  However, since this
    could cause confusion, the default is to make a copy and if a user
    wants a view, they must pass in an explicit flag.

    A VIEW of a vector uses the same underlying memory.
    WARNING: A View can be dangerous since multiple objects can access
    the same memory locations.
   */
  Vector(const Teuchos::RCP<Epetra_Vector>& source,
     NOX::Epetra::Vector::MemoryType memoryType =
     NOX::Epetra::Vector::CreateCopy,
     NOX::CopyType type = NOX::DeepCopy,
     Teuchos::RCP<NOX::Epetra::VectorSpace> vs =
     Teuchos::null);

  //! Construct by copying map and/or elements of an Epetra_Vector.
  /*! Allocates an entirely new vector. Does NOT allow for a view. */
  Vector(const Epetra_Vector& source,
     NOX::CopyType type = NOX::DeepCopy,
     Teuchos::RCP<NOX::Epetra::VectorSpace> vs =
     Teuchos::null);

  //! Copy constructor.
  Vector(const NOX::Epetra::Vector& source,
     NOX::CopyType type = NOX::DeepCopy);

  //! Destruct Vector.
  ~Vector();

  //@{ \name Access to underlying Petra vector.

  //! Get reference to underlying Epetra vector.
  virtual Epetra_Vector& getEpetraVector();

  //! Get const reference to underlying Epetra vector.
  virtual const Epetra_Vector& getEpetraVector() const;

  //@}

  //@{ \name Initialization methods.

  // derived
  virtual NOX::Abstract::Vector& init(double gamma);

  // derived
  virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed = 1);

  //! Copies source vector into "this".
  /*! NOTE: this will NOT copy the underlying vector space into the
    new vector. */
  virtual NOX::Abstract::Vector& operator=(const Epetra_Vector& y);

  // derived
  virtual NOX::Abstract::Vector& operator=(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

  // derived
  virtual NOX::Abstract::Vector& abs(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

  // derived
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

  //@}

  //@{ \name Update methods.

  // derived
  virtual NOX::Abstract::Vector& scale(double gamma);

  // derived
  virtual NOX::Abstract::Vector& scale(const NOX::Epetra::Vector& a);
  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a,
                 double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a,
                 double gamma = 0.0);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a,
                 double beta, const NOX::Epetra::Vector& b,
                 double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a,
                 double beta, const NOX::Abstract::Vector& b,
                 double gamma = 0.0);


  //@}

  //@{ \name Creating new Vectors.

  // derived
  virtual Teuchos::RCP<NOX::Abstract::Vector>
  clone(CopyType type = DeepCopy) const;

  /*!
   * \brief Create a MultiVector with \c numVecs+1 columns out of an array of
   * Vectors.  The vector stored under \c this will be the first column with
   * the remaining \c numVecs columns given by \c vecs.
   *
   * The implementation here creates a NOX::Epetra::MultiVector with
   * either Shape or Deep copies of the supplied vectors.
   */
  virtual Teuchos::RCP<NOX::Abstract::MultiVector>
  createMultiVector(const NOX::Abstract::Vector* const* vecs,
            int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

  /*!
   * \brief Create a MultiVector with \c numVecs columns.
   *
   * The implementation here creates a NOX::Epetra::MultiVector with
   * either Shape or Deep copies of the supplied vector.
   */
  virtual Teuchos::RCP<NOX::Abstract::MultiVector>
  createMultiVector(int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

  //@}

  //@{ \name Norms.

  // derived
  virtual double norm(NOX::Abstract::Vector::NormType type = TwoNorm) const;

  // derived
  virtual double norm(const NOX::Epetra::Vector& weights) const;
  virtual double norm(const NOX::Abstract::Vector& weights) const;

  //@}

  //@{ \name Inner products

  // derived
  virtual double innerProduct(const NOX::Epetra::Vector& y) const;
  virtual double innerProduct(const NOX::Abstract::Vector& y) const;

  //@}

  // derived
  virtual NOX::size_type length() const;

  // derived
  virtual void print(std::ostream& stream) const;

  //! Returns the NOX::Epetra::VectorSpace associated with this vector.
  virtual Teuchos::RCP<NOX::Epetra::VectorSpace>
  getVectorSpace() const;

 protected:

  //! Pointer to petra vector owned by this object
  Teuchos::RCP<Epetra_Vector> epetraVec;

  //! Pointer to the vector space.
  Teuchos::RCP<NOX::Epetra::VectorSpace> vectorSpace;

};
} // namespace Epetra
} // namespace NOX

#endif
