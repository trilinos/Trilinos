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

#ifndef LOCA_EXTENDED_MULTIVECTOR_H
#define LOCA_EXTENDED_MULTIVECTOR_H

#include "NOX_Common.H"  // for <iostream>, <vector>, and <string>
#include "NOX_Abstract_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Extended {
    class Vector;
  }
}

namespace LOCA {

  namespace Extended {

    /*! 
     * \brief Implemenatation of the NOX::Abstract::MultiVector class for 
     * extended multi-vectors comprised of an arbitrary number of 
     * multi-vectors and scalars.
     */
    /*
     * The implementations of the NOX::Abstract::MultiVector methods are 
     * defined in terms of the implementations of each stored abstract 
     * multi-vector and scalar.
     */

    class MultiVector : public NOX::Abstract::MultiVector {

      //! Declare LOCA::Extended::Vector as a friend class
      friend class Vector;

    public:

      //! Copy constructor.
      MultiVector(const MultiVector& source, 
		  NOX::CopyType type = NOX::DeepCopy);

      /*!
       * \brief Copy constructor that creates s sub copy or view of the 
       * given multivector
       */
      MultiVector(const MultiVector& source, 
		  vector<int>& index, bool view = false);

      //! Vector destructor.
      virtual ~MultiVector();

      //@{ \name Initialization methods.

      //! Initialize every element of this multi-vector with \c gamma.
      virtual NOX::Abstract::MultiVector& 
      init(double gamma);

      //! Initialize each element of this multi-vector with a random value
      virtual NOX::Abstract::MultiVector& 
      random(bool useSeed = false, int seed = 1);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector& 
      operator=(const NOX::Abstract::MultiVector& source);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual MultiVector& operator=(const MultiVector& y);

      /*! 
       * \brief Copy the vectors in \c source to a set of vectors in \c *this.
       * The \c index.size() vectors in \c source are copied to a subset of 
       * vectors in \c *this indicated by the indices given in \c index.
       */
      virtual NOX::Abstract::MultiVector&
      setBlock(const NOX::Abstract::MultiVector& source, vector<int>& index);
      virtual NOX::Abstract::MultiVector&
      setBlock(const MultiVector& source, vector<int>& index);

      //! Append the vectors in \c source to \c *this
      virtual NOX::Abstract::MultiVector&
      augment(const NOX::Abstract::MultiVector& source);
      virtual NOX::Abstract::MultiVector&
      augment(const MultiVector& source);

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
      update(double alpha, const MultiVector& a, 
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
      update(double alpha, const MultiVector& a, 
	     double beta, const MultiVector& b,
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
      update(double alpha, const MultiVector& a, 
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
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Creates a new multi-vector with \c numvecs columns
      virtual NOX::Abstract::MultiVector* 
      clone(int numvecs) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns whose 
       * columns are copies of the columns of \c *this given by \c index.
       */
      virtual NOX::Abstract::MultiVector* 
      subCopy(vector<int>& index) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns that 
       * shares the columns of \c *this given by \c index.
       */
      virtual NOX::Abstract::MultiVector* 
      subView(vector<int>& index) const;

      //@}

      //@{ \name Norms.

      //! Norm.
      virtual void
      norm(vector<double>& result,
	   NOX::Abstract::Vector::NormType type = 
	   NOX::Abstract::Vector::TwoNorm) const;

      //@}

      //@{ \name Dot products.

      /*! 
       * \brief Computes the dot product between each column of \c *this 
       * and \c y, scaled by \c alpha.
       */
      virtual void dot(double alpha, const NOX::Abstract::MultiVector& y,
		       NOX::Abstract::MultiVector::DenseMatrix& b) const;
      virtual void dot(double alpha, const MultiVector& y,
		       NOX::Abstract::MultiVector::DenseMatrix& b) const;
  
      //@}

      //! Return the length of multi-vector.
      virtual int length() const;

      //! Return the number of vectors in the multi-vector.
      virtual int numVectors() const;

      //! Print the vector.  This is meant for debugging purposes only 
      virtual void print() const;

      //! Sets the ith multi-vector
      virtual void setMultiVector(int i, const NOX::Abstract::MultiVector& v);

      //! Sets the scalar for row i, column j
      virtual void setScalar(int i, int j, double s);

      //! Returns const reference to the ith multi-vector
      virtual const NOX::Abstract::MultiVector& getMultiVector(int i) const;

      //! Returns reference to the ith multi-vector
      virtual NOX::Abstract::MultiVector& getMultiVector(int i);

      //! Returns const reference to the scalar for row i, column j
      virtual const double& getScalar(int i, int j) const;

      //! Returns reference to the scalar for row i, column j
      virtual double& getScalar(int i, int j);

    protected:

      /*!
       * \brief  Constructor that creates an empty extended multivector
       * holding nVecs multivectors, each multivector stores nColumns,
       * and there are nScalarRows rows of scalars.  Set view to true
       * when creating a multivector view.
       */
      MultiVector(int nColumns, int nVecs, int nScalarRows, bool view = false);

      //! Sets the pointer to the ith multivector
      void setMultiVectorPtr(int i, NOX::Abstract::MultiVector* v);

      //! Checks multi-vec argument dimensions are consistent
      void checkDimensions(const LOCA::Extended::MultiVector& a) const;

      //! Checks multi-vec argument dimensions are consistent
      void checkDimensions(const LOCA::Extended::MultiVector& a,
			   vector<int>& index) const;

      //! Checks multi-vec argument dimensions are consistent
      void checkDimensions(const LOCA::Extended::MultiVector& a,
		     const NOX::Abstract::MultiVector::DenseMatrix& b) const;

      //! Checks multi-vec argument dimensions are consistent for augment
      void checkAugmentDimensions(const LOCA::Extended::MultiVector& a) const;

      //! Checks validity of column index
      void checkIndex(int i) const;

      //! Checks validity of column and row index for scalars
      void checkIndex(int i, int j) const;

    protected:

      /*! 
       * \brief Number of columns in each multivec and number of scalar 
       * vector columns
       */
      int numColumns;

      //! Number of multivec block rows
      int numMultiVecRows;

      //! Number of scalar rows
      int numScalarRows;

      //! Array of multi-vector pointers, one for each block ROW
      vector< NOX::Abstract::MultiVector* > multiVectorPtrs;

      //! Array of scalars, one for each COLUMN
      /*!
       * We store the scalar vectors column wise instead of row wise
       * because for all but the matrix multiplication method, it is
       * more efficient.  This is also the only way to get sub-vector
       * views efficiently.
       */
      vector< NOX::Abstract::MultiVector::DenseMatrix* > scalarVectorPtrs;

      //! Pointers to each column as a LOCA::Extended::Vector
      mutable vector< LOCA::Extended::Vector* > extendedVectorPtrs;

      //! Flag indicating whether this multivec was created as a view
      bool isView;

    }; // class MultiVector
  } // namespace Extended
} // namespace LOCA

#endif
