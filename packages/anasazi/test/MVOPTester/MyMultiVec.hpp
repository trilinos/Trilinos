// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef MY_MULTIVECTOR_HPP
#define MY_MULTIVECTOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVec.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//! Simple example of a user's defined Anasazi::MultiVec class.
/*! 
 * This is a simple, single processor example of user's defined
 * MultiVec-derived class. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", or
 * "complex<double>".
 *
 * \author Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB) 
 *
 * \date Last modified on 01-Nov-05
 */
template <class ScalarType>
class MyMultiVec : public Anasazi::MultiVec<ScalarType>
{
public:

  //! Constructor for a \c NumberVecs vectors of length \c Length.
  MyMultiVec(const int Length, const int NumberVecs) :
    Length_(Length),
    NumberVecs_(NumberVecs)
  {
    Check();

    data_.resize(NumberVecs);
    ownership_.resize(NumberVecs);
    
    // Allocates memory to store the vectors.
    for (int v = 0 ; v < NumberVecs_ ; ++v)
    {
      data_[v] = new ScalarType[Length];
      ownership_[v] = true;
    }

    // Initializes all elements to zero.
    MvInit(0.0);
  }
  
  //! Constructor with already allocated memory
  MyMultiVec(const int Length, const std::vector<ScalarType*>& rhs) :
    Length_(Length),
    NumberVecs_(rhs.size())
  {
    Check();

    data_.resize(NumberVecs_);
    ownership_.resize(NumberVecs_);
    
    // Copies pointers from input array, set ownership so that
    // this memory is not free'd in the destructor
    for (int v = 0 ; v < NumberVecs_ ; ++v)
    {
      data_[v] = rhs[v];
      ownership_[v] = false;
    }
  }
  
  //! Copy constructor, performs a deep copy.
  MyMultiVec(const MyMultiVec& rhs) :
    Length_(rhs.GetGlobalLength()),
    NumberVecs_(rhs.NumberVecs_)
  {
    Check();

    data_.resize(NumberVecs_);
    ownership_.resize(NumberVecs_);
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
    {
      data_[v] = new ScalarType[Length_];
      ownership_[v] = true;
    }
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
    {
      for (int i = 0 ; i < Length_ ; ++i)
        (*this)(i, v) = rhs(i, v);
    }
  }
  
  //! Destructor
  ~MyMultiVec()
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      if (ownership_[v]) 
        delete[] data_[v];
  }
  
  //! Returns a clone of the current vector.
  MyMultiVec* Clone(const int NumberVecs) const
  {
    // FIXME
    MyMultiVec* tmp = new MyMultiVec(Length_, NumberVecs);
    
    //   for (int v = 0 ; v < NumberVecs ; ++v)
    //         for (int i = 0 ; i < Length_ ; ++i)
    //           (*tmp)(i, v) = (*this)(i, v);
    
    return(tmp);
  }
  
  // Returns a clone of the corrent multi-vector.
  MyMultiVec* CloneCopy() const
  {
    return(new MyMultiVec(*this));
  }
  
  //! Returns a clone copy of specified vectors.
  MyMultiVec* CloneCopy(const std::vector< int > &index) const
  {
    int size = index.size();
    MyMultiVec* tmp = new MyMultiVec(Length_, size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*tmp)(i, v) = (*this)(i, index[v]);
      }
    }
    
    return(tmp);
  }
  
  //! Returns a view of current vector (shallow copy)
  MyMultiVec* CloneViewNonConst(const std::vector< int > &index) 
  {
    int size = index.size();
    std::vector<ScalarType*> values(size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      values[v] = data_[index[v]];
    
    return(new MyMultiVec(Length_, values));
  }
  
  //! Returns a view of current vector (shallow copy), const version.
  const MyMultiVec* CloneView(const std::vector< int > &index) const
  {
    int size = index.size();
    std::vector<ScalarType*> values(size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      values[v] = data_[index[v]];
    
    return(new MyMultiVec(Length_, values));
  }
  
  ptrdiff_t GetGlobalLength () const
  {
    return(Length_);
  }
  
  int GetNumberVecs () const
  {
    return(NumberVecs_);
  }
  
  // Update *this with alpha * A * B + beta * (*this). 
  void MvTimesMatAddMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType> &A, 
                        const Teuchos::SerialDenseMatrix<int, ScalarType> &B, 
                        ScalarType beta)
  {
    
    assert (Length_ == A.GetGlobalLength());
    assert (B.numRows() == A.GetNumberVecs());
    assert (B.numCols() <= NumberVecs_);

    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert(MyA!=NULL);

    if ((*this)[0] == (*MyA)[0]) {
      // If this == A, then need additional storage ...
      // This situation is a bit peculiar but it may be required by
      // certain algorithms.
      
      std::vector<ScalarType> tmp(NumberVecs_);

      for (int i = 0 ; i < Length_ ; ++i) {
        for (int v = 0; v < A.GetNumberVecs() ; ++v) {
          tmp[v] = (*MyA)(i, v);
        }

        for (int v = 0 ; v < B.numCols() ; ++v) {
          (*this)(i, v) *= beta; 
          ScalarType res = Teuchos::ScalarTraits<ScalarType>::zero();

          for (int j = 0 ; j < A.GetNumberVecs() ; ++j) {
            res +=  tmp[j] * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
    else {
      for (int i = 0 ; i < Length_ ; ++i) {
        for (int v = 0 ; v < B.numCols() ; ++v) {
          (*this)(i, v) *= beta; 
          ScalarType res = 0.0;
          for (int j = 0 ; j < A.GetNumberVecs() ; ++j) {
            res +=  (*MyA)(i, j) * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
  }
  
  // Replace *this with alpha * A + beta * B. 
  void MvAddMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
                ScalarType beta,  const Anasazi::MultiVec<ScalarType>& B)
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    MyMultiVec* MyB;
    MyB = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(B)); 
    assert (MyB != 0);
    
    assert (NumberVecs_ == A.GetNumberVecs());
    assert (NumberVecs_ == B.GetNumberVecs());
    
    assert (Length_ == A.GetGlobalLength());
    assert (Length_ == B.GetGlobalLength());
    
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, v) = alpha * (*MyA)(i, v) + beta * (*MyB)(i, v);
      }
    }
  }
  
  // Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this). 
  void MvTransMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
                  Teuchos::SerialDenseMatrix< int, ScalarType >& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                  , Anasazi::ConjType conj
#endif
                 ) const
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetGlobalLength() == Length_);
    assert (NumberVecs_ <= B.numCols());
    assert (A.GetNumberVecs() <= B.numRows());
    
#ifdef HAVE_ANASAZI_EXPERIMENTAL
    if (conj == Anasazi::CONJ) {
#endif
      for (int v = 0 ; v < A.GetNumberVecs() ; ++v) {
        for (int w = 0 ; w < NumberVecs_ ; ++w) {
          ScalarType value = 0.0;
          for (int i = 0 ; i < Length_ ; ++i) {
            value += Teuchos::ScalarTraits<ScalarType>::conjugate((*MyA)(i, v)) * (*this)(i, w);
          }
          B(v, w) = alpha * value;
        }
      }
#ifdef HAVE_ANASAZI_EXPERIMENTAL
    } else {
      for (int v = 0 ; v < A.GetNumberVecs() ; ++v) {
        for (int w = 0 ; w < NumberVecs_ ; ++w) {
          ScalarType value = 0.0;
          for (int i = 0 ; i < Length_ ; ++i) {
            value += (*MyA)(i, v) * (*this)(i, w);
          }
          B(v, w) = alpha * value;
        }
      }
    }
#endif
  }
  
  
  // Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
  void MvDot (const Anasazi::MultiVec<ScalarType>& A, std::vector<ScalarType> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
              , Anasazi::ConjType conj
#endif
             ) const
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (NumberVecs_ <= (int)b.size());
    assert (NumberVecs_ == A.GetNumberVecs());
    assert (Length_ == A.GetGlobalLength());
    
#ifdef HAVE_ANASAZI_EXPERIMENTAL
    if (conj == Anasazi::CONJ) {
#endif
      for (int v = 0 ; v < NumberVecs_ ; ++v) {
        ScalarType value = 0.0;
        for (int i = 0 ; i < Length_ ; ++i) {
          value += (*this)(i, v) * Teuchos::ScalarTraits<ScalarType>::conjugate((*MyA)(i, v));
        }
        b[v] = value;
      }
#ifdef HAVE_ANASAZI_EXPERIMENTAL
    } else {
      for (int v = 0 ; v < NumberVecs_ ; ++v) {
        ScalarType value = 0.0;
        for (int i = 0 ; i < Length_ ; ++i) {
          value += (*this)(i, v) * (*MyA)(i, v);
        }
        b[v] = value;
      }
    }
#endif
  }  
  
	void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec ) const 
  {
    assert (NumberVecs_ <= (int)normvec.size());

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      MagnitudeType value = Teuchos::ScalarTraits<MagnitudeType>::zero();
      for (int i = 0 ; i < Length_ ; ++i) {
        MagnitudeType val = Teuchos::ScalarTraits<ScalarType>::magnitude((*this)(i, v));
        value += val * val;
      }
      normvec[v] = Teuchos::ScalarTraits<MagnitudeType>::squareroot(value);
    }
  }
  
  // Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
  // A are copied to a subset of vectors in *this indicated by the indices given 
  // in index.
  // FIXME: not so clear what the size of A and index.size() are...
  void SetBlock (const Anasazi::MultiVec<ScalarType>& A, 
                 const std::vector<int> &index)
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetNumberVecs() >= (int)index.size());
    assert (A.GetGlobalLength() == Length_);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, index[v])  = (*MyA)(i, v);
      }
    }
  }
  
  // Scale the vectors by alpha
  void MvScale( ScalarType alpha )
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, v) *= alpha;
      }
    }
  }

  // Scale the i-th vector by alpha[i]
  void MvScale( const std::vector<ScalarType>& alpha )
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, v) *= alpha[v];
      }
    }
  }

  // Fill the vectors in *this with random numbers.
  void MvRandom ()
  {
    Teuchos::SerialDenseMatrix<int,ScalarType> R( Length_, NumberVecs_ );
    Teuchos::randomSyncedMatrix( R );
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, v) = R(i, v);
      }
    }
  }

  // Replace each element of the vectors in *this with alpha.
  void  MvInit (ScalarType alpha)
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v) {
      for (int i = 0 ; i < Length_ ; ++i) {
        (*this)(i, v) = alpha;
      }
    }
  }
  
  void MvPrint (std::ostream &os) const
  {
    os << "Object MyMultiVec" << std::endl;
    os << "Number of rows = " << Length_ << std::endl;
    os << "Number of vecs = " << NumberVecs_ << std::endl;
    
    for (int i = 0 ; i < Length_ ; ++i)
      {
        for (int v = 0 ; v < NumberVecs_ ; ++v)
          os << (*this)(i, v) << " ";
        os << std::endl;
      }
  }
  
  inline ScalarType& operator()(const int i, const int j)
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    // 
    return(data_[j][i]);
  }
  
  inline const ScalarType& operator()(const int i, const int j) const
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    // 
    return(data_[j][i]);
  }
  
  ScalarType* operator[](int v)
  {
    if (v < 0 || v >= NumberVecs_) throw(-1);
    return(data_[v]);
  }
  
  ScalarType* operator[](int v) const
  {
    return(data_[v]);
  }
  
private:
  void Check()
  {
    if (Length_ <= 0)
      throw("Length must be positive");

    if (NumberVecs_ <= 0)
      throw("Number of vectors must be positive");
  }

  //! Length of the vectors
  const int Length_;
  //! Number of multi-vectors
  const int NumberVecs_;
  //! Pointers to the storage of the vectors.
  std::vector<ScalarType*> data_;
  //! If \c true, then this object owns the vectors and must free them in dtor.
  std::vector<bool> ownership_;

};

/*! \relates MyMultiVec
    Output stream operator for handling the printing of MyMultiVec.
*/
template <typename ScalarType>
std::ostream& operator<<(std::ostream& os, const MyMultiVec<ScalarType>& Obj)
{
  Obj.MvPrint(os);
  return os;
}

#endif // MY_MULTIVECTOR_HPP
