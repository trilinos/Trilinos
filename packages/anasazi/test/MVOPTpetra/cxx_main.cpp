//@HEADER
// ************************************************************************
// 
//
//                 Anasazi: Block Eigensolvers Package
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
// ************************************************************************
//@HEADER
//
//  This test instantiates the Anasazi classes using a complex scalar type
//  and checks functionality.
//

#include "AnasaziConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

typedef int OrdinalType;
#if 0
#ifdef HAVE_COMPLEX
typedef std::complex<double> ScalarType;
#elif HAVE_COMPLEX_H
typedef ::complex<double> ScalarType;
#endif
#endif
typedef double ScalarType;

// Gets the type of magnitude type (for example, the magnitude type
// of "complex<double>" is "double").
#include "AnasaziOperator.hpp"
#include "AnasaziMultiVec.hpp"

// ==============================================================
// Simple class to extend Tpetra::Vector into Tpetra::MultiVector
//
// This class stores a std::vector of RefCountPtr's to allocated
// Tpetra::Vector's. It overloads the (i, j) operator so that
// one can easily access the j-th element of the i-th vector.
// It also defines some other basic methods, mainly to fulfill
// what Anasazi wants.
// ==============================================================

namespace Tpetra {

  template<class OrdinalType, class ScalarType>
  class MultiVector {

    public:

      // Basic constructor
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, const OrdinalType NumVectors) :
        vectorSpace_(vectorSpace),
        NumVectors_(NumVectors)
      {
        Init(); 

        array_.resize(getNumVectors());

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (vectorSpace));
        }
      }

      // Creates a deep copy of the input set of vectors.
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                  std::vector<Tpetra::Vector<OrdinalType, ScalarType> const *> list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init(); 

        array_.resize(getNumVectors());

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // deep copy of each of the vectors
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (*(list[i])));
        }
      }

      // Creates a shallow copy of the input set of vectors.
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                  std::vector<Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > > list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init();

        array_.resize(NumVectors_);

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // copy RefCountPtr's from the list to this array.
          array_[i] = list[i];
        }
      }

      // Copy constructor.
      MultiVector(const MultiVector& rhs) :
        vectorSpace_(rhs.vectorSpace()),
        NumVectors_(rhs.getNumVectors())
      {
        Init();

        array_.resize(NumVectors_);

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (*(rhs.GetVector(i))));
        }
      }

      //! Returns the global number of entries.
      OrdinalType getNumGlobalEntries() const
      {
        return(vectorSpace_.getNumGlobalEntries());
      }

      //! Returns the number of entries on the calling image.
      OrdinalType getNumMyEntries() const
      {
        return(vectorSpace_.getNumMyEntries());
      }

      //! Returns the number of vectors in this multivector.
      OrdinalType getNumVectors() const
      {
        return(NumVectors_);
      }

      //! Returns a reference to the vector space of this multivector.
      VectorSpace<OrdinalType, ScalarType> const& vectorSpace () const
      {
        return(vectorSpace_);
      }

      //! Returns a reference to the i-th element of the j-th vector.
      ScalarType& operator() (const OrdinalType i, const OrdinalType j)
      {
        return((*array_[j])[i]);
      }

      //! Returns a reference to the i-th element of the j-th vector (const version)
      ScalarType const& operator() (const OrdinalType i, const OrdinalType j) const
      {
        return((*array_[j])[i]);
      }

      //! Sets all elements of all vectors to the given value.
      void setAllToScalar(ScalarType const value)
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          array_[i]->setAllToScalar(value);
      }

      //! Sets all elements of all vectors to random value.
      void setAllToRandom()
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          array_[i]->setAllToRandom();
      }

      //! Prints the vector to cout. FIXME
      void Print() const
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          cout << (*array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector.
      Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > GetRCP(const int i)
      {
        return(array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector (const version).
      Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> const > GetRCP(const int i) const
      {
        return(array_[i]);
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector.
      Tpetra::Vector<OrdinalType, ScalarType>* GetVector(const int i)
      {
        return(array_[i].get());
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector (const version).
      Tpetra::Vector<OrdinalType, ScalarType> const* GetVector(const int i) const
      {
        return(array_[i].get());
      }

      void norm1(ScalarType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = array_[i]->norm1();
        }
      }

      void norm2(typename Teuchos::ScalarTraits<ScalarType>::magnitudeType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = ScalarOne_ * (array_[i]->norm2()); // FIXME: DOES NOT WORK WITH COMPLEX
        }
      }

      void normInf(ScalarType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = array_[i]->normInf();
        }
      }

    private:

      void Init()
      {
        OrdinalZero_ = Teuchos::ScalarTraits<OrdinalType>::zero();
        OrdinalOne_  = Teuchos::ScalarTraits<OrdinalType>::one();

        ScalarZero_ = Teuchos::ScalarTraits<ScalarType>::zero();
        ScalarOne_  = Teuchos::ScalarTraits<ScalarType>::one();
      }

      std::vector<Teuchos::RefCountPtr<Vector<OrdinalType, ScalarType> > > array_;
      const VectorSpace<OrdinalType, ScalarType>& vectorSpace_;
      OrdinalType NumVectors_;

      OrdinalType OrdinalZero_;
      OrdinalType OrdinalOne_;

      ScalarType ScalarZero_;
      ScalarType ScalarOne_;

  }; // class MultiVector

} // namespace Tpetra

namespace Anasazi {

template<class OrdinalType, class ScalarType>
class TpetraOperator : public Anasazi::Operator<ScalarType>
{
public:
  TpetraOperator(Tpetra::CisMatrix<OrdinalType, ScalarType>& Op) :
    Op_(Op)
  {}

  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
			    Anasazi::MultiVec<ScalarType>& Y) const
  {
    Tpetra::MultiVector<OrdinalType, ScalarType>* MyX;
    MyX = dynamic_cast<Tpetra::MultiVector<OrdinalType, ScalarType>*>(const_cast<Anasazi::MultiVec<ScalarType>*>(&X)); 
    assert (MyX != 0);
    
    Tpetra::MultiVector<OrdinalType, ScalarType>* MyY;
    MyY = dynamic_cast<Tpetra::MultiVector<OrdinalType, ScalarType>*>(&Y); 
    assert (MyY != 0);

    for (int i = 0 ; i < X.GetNumberVecs() ; ++i)
      Op_.apply(*(MyX->GetVector(i)), *(MyY->GetVector(i)), false);

    return(Anasazi::Ok);
  }

private:
  Tpetra::CisMatrix<OrdinalType, ScalarType>& Op_;

}; // class TpetraOperator

class TpetraMultiVec : public Anasazi::MultiVec<ScalarType>, 
                       public Tpetra::MultiVector<OrdinalType, ScalarType>
{
public:
  //@{ \name Constructors and Destructors
  //! Basic constructor
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, const int numvecs) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, numvecs)
  {}

  // This is a deep copy
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                 std::vector<Tpetra::Vector<OrdinalType, ScalarType> const *> list) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, list)
  {}

  // This is a shalloe copy
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                 std::vector<Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > > list) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, list)
  {}

  //! Copy constructor.
  TpetraMultiVec(const Tpetra::MultiVector<OrdinalType, ScalarType>& rhs) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(rhs)
  {}

  //! Destructor
  virtual ~TpetraMultiVec() {};

  //@}
  //@{ \name Creation methods

  /*! \brief Creates a new empty EpetraMultiVec containing \c numvecs columns.

    \returns Pointer to an EpetraMultiVec
  */
  MultiVec<ScalarType> * Clone (const int numvecs) const
  {
    return(new TpetraMultiVec(vectorSpace(), numvecs));
  }

  /*! \brief Creates a new EpetraMultiVec and copies contents of \c *this into
      the new vector (deep copy).

    \returns Pointer to an EpetraMultiVec
  */	
  MultiVec<ScalarType> * CloneCopy () const
  {
    return(new TpetraMultiVec(*this));
  }

  /*! \brief Creates a new EpetraMultiVec and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
    The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.

    \returns Pointer to an EpetraMultiVec
  */
  MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const
  {
    std::vector<Tpetra::Vector<OrdinalType, ScalarType> const*> list(index.size());
    for (int i = 0 ; i < index.size() ; ++i)
      list[i] = this->GetVector(i);

    return(new TpetraMultiVec(vectorSpace(), list));
  }
    
  /*! \brief Creates a new EpetraMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraMultiVec
  */
  MultiVec<ScalarType> * CloneView ( const std::vector<int>& index )
  {
    std::vector<Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > > list(index.size());
    for (int i = 0 ; i < index.size() ; ++i)
      list[i] = this->GetRCP(i);

    return(new TpetraMultiVec(vectorSpace(), list));
  }

  //@}
  //@{ \name Attribute methods	

  //! Obtain the vector length of *this.
  int GetNumberVecs () const { return getNumVectors(); }

  //! Obtain the number of vectors in *this.
  int GetVecLength () const { return getNumGlobalEntries(); }

  //@}

  //@{ \name Update methods
  /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
   */
  void MvTimesMatAddMv(const ScalarType alpha, const MultiVec<ScalarType>& A, 
                       const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta )
  {
    // FIXME: only works in serial
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert(MyA!=NULL);

    if (this == MyA)
    {
      // If this == A, then need additional storage ...
      // This situation is a bit peculiar but it may be required by
      // certain algorithms.
      
      std::vector<ScalarType> tmp(getNumVectors());

      for (int i = 0 ; i < getNumMyEntries() ; ++i)
      {
	for (int v = 0; v < getNumMyEntries() ; ++v) tmp[v] = (*MyA)(i, v);

        for (int v = 0 ; v < B.numCols() ; ++v)
        {
          (*this)(i, v) *= beta; 
          ScalarType res = 0.0;
          for (int j = 0 ; j < A.GetNumberVecs() ; ++j)
          {
            res +=  tmp[j] * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
    else
    {
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
      {
        for (int v = 0 ; v < B.numCols() ; ++v)
        {
          (*this)(i, v) *= beta; 
          ScalarType res = 0.0;
          for (int j = 0 ; j < A.GetNumberVecs() ; ++j)
          {
            res +=  (*MyA)(i, j) * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
  }

  /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
   */
  void MvAddMv (const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta,
                const MultiVec<ScalarType>& B)
  {
    // FIXME: only works in serial
    
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    TpetraMultiVec* MyB;
    MyB = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(B)); 
    assert (MyB != 0);
    
    for (int v = 0 ; v < getNumVectors() ; ++v)
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
	(*this)(i, v) = alpha * (*MyA)(i, v) + beta * (*MyB)(i, v);
  }

  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
  */
  void MvTransMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B ) const
  {
    // FIXME: only works in serial

    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetVecLength() == getNumMyEntries());
    assert (getNumVectors() == B.numCols());
    assert (A.GetNumberVecs() == B.numRows());
    
    for (int v = 0 ; v < A.GetNumberVecs() ; ++v)
    {
      for (int w = 0 ; w < getNumVectors() ; ++w)
      {
        ScalarType value = 0.0;
        for (int i = 0 ; i < getNumMyEntries() ; ++i)
        {
          value += (*MyA)(i, v) * (*this)(i, w);
        }
        B(v, w) = alpha * value;
      }
    }
  }

  /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^T(this[i])\f$ where \c A[i] is the i-th column of \c A.
  */
  void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType>* b ) const
  {
    // FIXME: only works in serial

    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (getNumVectors() == (int)b->size());
    assert (getNumVectors() == A.GetNumberVecs());
    assert (getNumMyEntries() == A.GetVecLength());
    
    for (int v = 0 ; v < getNumVectors() ; ++v)
    {
      ScalarType value = 0.0;
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
        value += (*this)(i, v) * (*MyA)(i, v);
      (*b)[v] = value;
    }
  }

  //@}
  //@{ \name Norm method

  /*! \brief Compute the 2-norm of each individual vector of \c *this.  
    Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
  void MvNorm ( std::vector<Teuchos::ScalarTraits<ScalarType>::magnitudeType >* normvec ) const 
  {
    assert (normvec != 0);
    assert (getNumVectors() == (int)normvec->size());
    
    norm2(&((*normvec)[0]));
  }
  //@}

  //@{ \name Initialization methods
  /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
  */
  void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index )
  {
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetNumberVecs() >= (int)index.size());
    assert (A.GetVecLength() == getNumMyEntries());
    
    cout << getNumMyEntries()  << endl;

    for (unsigned int v = 0 ; v < index.size() ; ++v)
    {
      cout << index[v] << " .. " << getNumVectors() << endl;
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
	(*this)(i, index[v])  = (*MyA)(i, v);
    }
  }

  /*! \brief Fill the vectors in \c *this with random numbers.
   */
  void MvRandom()
  {
    setAllToRandom();
  }

  /*! \brief Replace each element of the vectors in \c *this with \c alpha.
   */
  void MvInit (const ScalarType alpha)
  {
    setAllToScalar(alpha);
  }

  //@}
  //@{ \name Print method.
  /*! \brief Print \c *this EpetraMultiVec.
   */
  void MvPrint( ostream& os ) const 
  {
    Print(); // FIXME
  }
  //@}


}; // class TpetraMultiVec
}; // namespace Anasazi

// =========== //
// main driver //
// =========== //

#include "AnasaziBlockJacobiDavidson.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockJacobiDavidson.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 10;
  OrdinalType indexBase = OrdinalZero;

  // 1) Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // 2) We can now create a space:

  Tpetra::ElementSpace<OrdinalType> elementSpace(length, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> 
    vectorSpace(elementSpace, platformV);

  // 3) We now setup a matrix, which is diagonal.
  //    To that aim, we need to extract the list of locally owned
  //    ID's from elementSpace.
  
  OrdinalType NumMyElements = elementSpace.getNumMyElements();
  vector<OrdinalType> MyGlobalElements = elementSpace.getMyGlobalElements();
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(vectorSpace);

  // stencil on a given row is [b a c]
  ScalarType a = ScalarOne * 2;
  ScalarType b = ScalarOne * (-1);
  ScalarType c = ScalarOne * (-1);

  for (OrdinalType LID = OrdinalZero ; LID < NumMyElements ; ++LID)
  {
    OrdinalType GID = MyGlobalElements[LID];
    if (GID != OrdinalZero)
      matrix.submitEntry(Tpetra::Add, GID, b, GID - 1);
    if (GID != length - 1)
      matrix.submitEntry(Tpetra::Add, GID, c, GID + 1);
    matrix.submitEntry(Tpetra::Add, GID, a, GID);
  }

  matrix.fillComplete();

  //  ....

  Tpetra::CisMatrix<OrdinalType,ScalarType> matrixB(vectorSpace);

  for (OrdinalType LID = OrdinalZero ; LID < NumMyElements ; ++LID)
  {
    OrdinalType GID = MyGlobalElements[LID];
    matrixB.submitEntry(Tpetra::Add, GID, a, GID);
  }

  matrixB.fillComplete();

  OrdinalType NumVectors = OrdinalOne * 3;

  Teuchos::RefCountPtr<Anasazi::TpetraOperator<OrdinalType, ScalarType> > A = 
    Teuchos::rcp(new Anasazi::TpetraOperator<OrdinalType, ScalarType>(matrix));
  Teuchos::RefCountPtr<Anasazi::TpetraOperator<OrdinalType, ScalarType> > B = 
    Teuchos::rcp(new Anasazi::TpetraOperator<OrdinalType, ScalarType>(matrixB));
  Teuchos::RefCountPtr<Anasazi::TpetraOperator<OrdinalType, ScalarType> > K = 
    Teuchos::rcp(new Anasazi::TpetraOperator<OrdinalType, ScalarType>(matrix));






  Anasazi::ReturnType returnCode = Anasazi::Ok;	

  ScalarType     TOL    = 1e-8;

  int MyPID = 0;

  int MAXITER = 500;
  int NEV = 3; 
  int SMIN = 20;
  int SMAX = 30;
  // FIXME: if I set BLOCKSIZE = 2 then nothing good happens...
  int BLOCKSIZE = 1;

  ScalarType TARGET = 0.5;

  // ================== //
  // Sets up the solver //
  // ================== //

  typedef Anasazi::MultiVec<ScalarType> MV;        
  typedef Anasazi::Operator<ScalarType> OP;        
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  //typedef Anasazi::MultiVecTraits<ScalarType, MyMultiVec> MVT;

  Teuchos::RefCountPtr<MV> ivec = Teuchos::rcp(new Anasazi::TpetraMultiVec(vectorSpace, BLOCKSIZE));        
  ivec->MvRandom();

  // Create the eigenproblem.
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<ScalarType, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<ScalarType, MV, OP>(A, B, ivec) );

  //MyProblem->SetPrec(K); 

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true); 

  // Set the number of eigenvalues requested
  MyProblem->SetNEV(NEV);

  // Inform the eigenproblem that you are finishing passing it information
  MyProblem->SetProblem();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<ScalarType> > MyOM =
    Teuchos::rcp(new Anasazi::OutputManager<ScalarType>(MyPID));
  MyOM->SetVerbosity(Anasazi::FinalSummary);	

  // Create a sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<ScalarType, MV, OP> > MySM =
    Teuchos::rcp(new Anasazi::BasicSort<ScalarType, MV, OP>("SM"));


  // Create parameter list to pass into solver
  // FIXME: ADD PARAMTERS
  Teuchos::ParameterList MyPL;
  MyPL.set("Block Size", BLOCKSIZE);
  MyPL.set("SMIN", SMIN);
  MyPL.set("SMAX", SMAX);
  MyPL.set("Max Iters", MAXITER);
  MyPL.set("Tol", TOL);
  MyPL.set("Target", TARGET);

  // Initialize the Block Jacobi-Davidson solver
  Anasazi::BlockJacobiDavidson<ScalarType, MagnitudeType, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);
                           
  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  MySolver.currentStatus();

}
