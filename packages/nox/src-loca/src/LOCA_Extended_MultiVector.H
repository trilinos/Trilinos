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
       * \brief Copy constructor that creates a new multivector with
       * \em nColumns columns.
       */
      MultiVector(const MultiVector& source, int nColumns);

      /*!
       * \brief Copy constructor that creates a sub copy or view of the 
       * given multivector
       */
      MultiVector(const MultiVector& source, vector<int>& index, bool view);

      //! Vector destructor.
      virtual ~MultiVector();

      //@{ \name Implementation of NOX::Abstract::MultiVector methods.

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
       * \c a is a multivector, \c b is a dense matrix, \c x = \c *this,
       * and op(b) = b if transb = Teuchos::NO_TRANS and op(b) is b transpose
       * if transb = Teuchos::TRANS.
       */
      virtual NOX::Abstract::MultiVector& 
      update(Teuchos::ETransp transb, double alpha, 
	     const NOX::Abstract::MultiVector& a, 
	     const NOX::Abstract::MultiVector::DenseMatrix& b, 
	     double gamma = 0.0);
      virtual NOX::Abstract::MultiVector& 
      update(Teuchos::ETransp transb, double alpha, const MultiVector& a, 
	     const NOX::Abstract::MultiVector::DenseMatrix& b, 
	     double gamma = 0.0);

      /*! 
       * \brief Create a new %MultiVector of the same underlying type by
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

      //! Norm.
      virtual void
      norm(vector<double>& result,
	   NOX::Abstract::Vector::NormType type = 
	   NOX::Abstract::Vector::TwoNorm) const;

      /*! 
       * \brief Computes the matrix-matrix product \f$\alpha * y^T * (*this)\f$
       */
      virtual void multiply(double alpha, const NOX::Abstract::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const;
      virtual void multiply(double alpha, const MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const;

      //! Return the length of multi-vector.
      virtual int length() const;

      //! Return the number of vectors in the multi-vector.
      virtual int numVectors() const;

      //! Print the vector.  This is meant for debugging purposes only 
      virtual void print() const;

      //@}

      //@{ \name Additional public accessor methods

      //! Returns const reference to the ith multi-vector
      virtual const NOX::Abstract::MultiVector& getMultiVector(int i) const;

      //! Returns reference to the ith multi-vector
      virtual NOX::Abstract::MultiVector& getMultiVector(int i);

      //! Returns const reference to scalar matrix
      virtual const NOX::Abstract::MultiVector::DenseMatrix& 
      getScalars() const;

      //! Returns reference to scalar matrix
      virtual NOX::Abstract::MultiVector::DenseMatrix& getScalars();

      //! Returns const reference to the scalar for row i, column j
      virtual const double& getScalar(int i, int j) const;

      //! Returns reference to the scalar for row i, column j
      virtual double& getScalar(int i, int j);

      /*!
       * \brief Return a reference to the i-th column of the multivector as
       * an abstract vector.  
       */
      virtual LOCA::Extended::Vector& getVector(int i);
    
      /*!
       * \brief Return a const reference to the i-th column of the 
       * multivector as an abstract vector.  
       */
      virtual const LOCA::Extended::Vector& getVector(int i) const;

      //! Returns number of scalars rows
      virtual int getNumScalarRows() const;

      //! Returns number of multi vectors
      virtual int getNumMultiVectors() const;

      //@}

    protected:

      //! Constructor that creates an empty multivector to be filled in later.
      MultiVector(int nColumns, int nVectorRows, int nScalarRows);

      //@{ \name Additional interface derived classes must implement

      //! Generate a derived extended vector
      /*!
       * Extended multi-vectors derived from this class should implement
       * this method and return a vector of the appropriate type.
       * This allows the operator[] methods to work correctly for 
       * derived classes.
       */
      virtual LOCA::Extended::Vector* 
      generateVector(int nVecs, int nScalarRows) const;

      //@}

      //@{ \name Other protected helper functions

      //! Sets the pointer to the ith multivector
      void setMultiVectorPtr(int i, NOX::Abstract::MultiVector* v);

      //! Checks multi-vec argument dimensions are consistent
      void checkDimensions(const string& callingFunction, 
			   const LOCA::Extended::MultiVector& a) const;

      //! Checks validity of column index
      void checkIndex(const string& callingFunction, int i) const;

      //! Checks validity of column and row index for scalars
      void checkIndex(const string& callingFunction, int i, int j) const;

      //! Checks is index array is contiguous
      bool isContiguous(const vector<int>& index) const;

      //@}

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

      //! Dense matrix of scalars
      NOX::Abstract::MultiVector::DenseMatrix* scalarsPtr;

      //! Pointers to each column as a LOCA::Extended::Vector
      mutable vector< LOCA::Extended::Vector* > extendedVectorPtrs;

      //! Flag indicating whether this vector is a view
      bool isView;

    }; // class MultiVector
  } // namespace Extended
} // namespace LOCA

#endif
