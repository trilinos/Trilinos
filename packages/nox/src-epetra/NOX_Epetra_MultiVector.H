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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_MULTIVECTOR_H
#define NOX_EPETRA_MULTIVECTOR_H

#include "NOX_Abstract_MultiVector.H" // base class

#ifdef HAVE_NOX_MULTIVECS

// Forward declaration
class Epetra_MultiVector;

namespace NOX {
  namespace Epetra {
    class Vector;
  }
}

namespace NOX {

  namespace Epetra {

    //! Implementation of NOX::Abstract::MultiVector for %Epetra multi-vectors.
    class MultiVector : public virtual NOX::Abstract::MultiVector {

    public:			

      //! Construct by copying map and/or elements of an Epetra_MultiVector.
      /*! Option for creating a "View": Construct a NOX::Epetra::MultiVector 
	using the SAME memory location used in the source Epetra_MultiVector.
	WARNING: A View can be dangerous since multiple objects can access
	the same memory locations.  
      */
      MultiVector(Epetra_MultiVector& source, 
		  NOX::CopyType type = NOX::DeepCopy,
		  bool createView = false);

      //! Construct by copying map and/or elements of an Epetra_MultiVector.
      MultiVector(const Epetra_MultiVector& source, 
		  NOX::CopyType type = NOX::DeepCopy);

      //! Copy constructor.
      MultiVector(const NOX::Epetra::MultiVector& source, 
		  NOX::CopyType type = NOX::DeepCopy);

      //! Destruct MultiVector.
      ~MultiVector();

      //@{ \name Access to underlying Petra vector.

      //! Get reference to underlying Epetra vector.
      virtual Epetra_MultiVector& getEpetraMultiVector();

      //! Get const reference to underlying Epetra vector.
      virtual const Epetra_MultiVector& getEpetraMultiVector() const;
 
      //@}

      //@{ \name Initialization methods.

      //! Initialize every element of this multi-vector with \c gamma.
      virtual NOX::Abstract::MultiVector& 
      init(double value);

      //! Initialize each element of this multi-vector with a random value
      virtual NOX::Abstract::MultiVector& 
      random(bool useSeed = false, int seed = 1);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector& 
      operator=(const Epetra_MultiVector& source);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector& 
      operator=(const NOX::Epetra::MultiVector& source);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector& 
      operator=(const NOX::Abstract::MultiVector& source);

      /*! 
       * \brief Copy the vectors in \c source to a set of vectors in \c *this.
       * The \c index.size() vectors in \c source are copied to a subset of 
       * vectors in \c *this indicated by the indices given in \c index.
       */
      virtual NOX::Abstract::MultiVector&
      setBlock(const NOX::Abstract::MultiVector& source, 
	       vector<int>& index);
      virtual NOX::Abstract::MultiVector&
      setBlock(const NOX::Epetra::MultiVector& source, 
	      vector<int>& index);

      //! Append the vectors in \c source to \c *this
      virtual NOX::Abstract::MultiVector&
      augment(const NOX::Abstract::MultiVector& source);
      virtual NOX::Abstract::MultiVector&
      augment(const NOX::Epetra::MultiVector& source);

      //@}

      //@{ \name Accessor methods

      /*!
       * \brief Return a reference to the i-th column of the multivector as
       * an abstract vector.  
       */
      virtual NOX::Abstract::Vector& operator [] (int i);
    
      /*!
       * \brief Return a const reference to the i-th column of the 
       * multivector as an abstract vector.  
       */
      virtual const NOX::Abstract::Vector& operator [] (int i) const;

      //@}

      //@{ \name Update methods.

      /*! 
       * \brief Compute x = (alpha * a) + (gamma * x) where \c a is a 
       * multi-vector and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Abstract::MultiVector& a, 
	     double gamma = 0.0);
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Epetra::MultiVector& a, 
	     double gamma = 0.0);

      /*! 
       * \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where
       * \c a and \c b are multi-vectors and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Abstract::MultiVector& a, 
	     double beta, const NOX::Abstract::MultiVector& b,
	     double gamma = 0.0);
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Epetra::MultiVector& a, 
	     double beta, const NOX::Epetra::MultiVector& b,
	     double gamma = 0.0);

      /*! 
       * \brief Compute x = (alpha * a * b) + (gamma * x) where 
       * \c a is a multivector, \c b is a dense matrix, and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Abstract::MultiVector& a, 
	     const NOX::Abstract::MultiVector::DenseMatrix& b, 
	     double gamma = 0.0);
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Epetra::MultiVector& a, 
	     const NOX::Abstract::MultiVector::DenseMatrix& b, 
	     double gamma = 0.0);

      //@}

      //@{ \name Creating new MultiVectors. 

      /*! 
       * \brief Create a new %Vector of the same underlying type by
       * cloning "this", and return a pointer to the new vector.  
       *
       * If type is NOX::DeepCopy, then we need to create an exact replica
       * of "this". Otherwise, if type is NOX::ShapeCopy, we need only
       * replicate the shape of "this". Note that there is 
       * <em>no assumption</em> that a vector created by ShapeCopy is 
       * initialized to zeros.
       * 
       * \return Pointer to newly created vector or NULL if clone is not 
       * supported. 
       */
      virtual NOX::Abstract::MultiVector* 
      clone(CopyType type = DeepCopy) const;

      //! Creates a new multi-vector with \c numvecs columns
      virtual NOX::Abstract::MultiVector* clone(int numvecs) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns whose 
       * columns are copies of the columns of \c *this given by \c index.
       */
      virtual NOX::Abstract::MultiVector* 
      subCopy(vector<int>& index) const;

      /*! 
       * \brief Creates a new multi-vector with \c ndex.size() columns that 
       * shares the columns of \c *this given by \c index.
       */
      virtual NOX::Abstract::MultiVector* 
      subView(vector<int>& index) const;

      //@}

      //@{ \name Norms.

      //! Norm.
      virtual void
      norm(vector<double>& result, NOX::Abstract::Vector::NormType type = 
	   NOX::Abstract::Vector::TwoNorm) const;

      //@}

      //@{ \name Dot products

      /*! 
       * \brief Computes the dot product between each column of \c *this 
       * and \c y, scaled by \c alpha.
       */
      virtual void dot(double alpha, const NOX::Abstract::MultiVector& y,
		       NOX::Abstract::MultiVector::DenseMatrix& b) const;
      virtual void dot(double alpha, const NOX::Epetra::MultiVector& y,
		       NOX::Abstract::MultiVector::DenseMatrix& b) const;

      //@}

      //! Return the length of multi-vector.
      virtual int length() const;

      //! Return the number of vectors in the multi-vector.
      virtual int numVectors() const;

      //! Print the vector.  This is meant for debugging purposes only 
      virtual void print() const;

    protected:

       //! Pointer to petra vector owned by this object
      Epetra_MultiVector* epetraMultiVec;

      //! NOX::Epetra::Vector's for each column of the multivector
      /*!
       * Each Epetra_Vector in the NOX::Epetra::Vector has a view into a 
       * column of the multivector and get filled in as needed by operator[].
       */
      mutable vector<NOX::Epetra::Vector*> noxEpetraVectors;

      //! Constructor (Protected)
      MultiVector(int numvecs);

      //! Checks whether an index is valid.  Throws an error if invalid
      void checkIndex(int idx) const;

    };
  } // namespace Epetra
} // namespace NOX

#endif // HAVE_NOX_MULTIVECS

#endif // NOX_EPETRA_MULTIVECTOR_H
