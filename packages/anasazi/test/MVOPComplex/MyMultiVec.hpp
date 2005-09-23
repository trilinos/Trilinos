#ifndef MY_MULTIVECTOR_HPP
#define MY_MULTIVECTOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"


class MyMultiVec : public Anasazi::MultiVec<ScalarType>
{

  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType magnitudeType;

  public:

  MyMultiVec(const int Length, const int NumberVecs) :
    Length_(Length),
    NumberVecs_(NumberVecs)
  {
    data_.resize(NumberVecs);
    ownership_.resize(NumberVecs);
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      {
        data_[v] = new ScalarType[Length];
        ownership_[v] = true;
      }
    
    MvInit(0.0);
  }
  
  MyMultiVec(const int Length, const std::vector<ScalarType*>& rhs) :
    Length_(Length),
    NumberVecs_(rhs.size())
  {
    data_.resize(NumberVecs_);
    ownership_.resize(NumberVecs_);
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      {
        data_[v] = rhs[v];
        ownership_[v] = false;
      }
  }
  
  MyMultiVec(const MyMultiVec& rhs) :
    Length_(rhs.GetVecLength()),
    NumberVecs_(rhs.NumberVecs_)
  {
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
  
  ~MyMultiVec()
  {
    // FIXME!
    //        for (int v = 0 ; v < NumberVecs_ ; ++v)
    //        if (ownership_[v]) delete[] data_[v];
  }
  
  MyMultiVec* Clone(const int NumberVecs) const
  {
    MyMultiVec* tmp = new MyMultiVec(Length_, NumberVecs);
    
    //   for (int v = 0 ; v < NumberVecs ; ++v)
    //         for (int i = 0 ; i < Length_ ; ++i)
    //           (*tmp)(i, v) = (*this)(i, v);
    
    return(tmp);
  }
  
  MyMultiVec* CloneCopy() const
  {
    return(new MyMultiVec(*this));
  }
  
  MyMultiVec* CloneCopy(const std::vector< int > &index) const
  {
    int size = index.size();
    MyMultiVec* tmp = new MyMultiVec(Length_, size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      for (int i = 0 ; i < Length_ ; ++i)
	(*tmp)(i, v) = (*this)(i, index[v]);
    
    return(tmp);
  }
  
  MyMultiVec* CloneView(const std::vector< int > &index) 
  {
    int size = index.size();
    std::vector<ScalarType*> values(size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      values[v] = data_[index[v]];
    
    return(new MyMultiVec(Length_, values));
  }
  
  const MyMultiVec* CloneView(const std::vector< int > &index) const
  {
    int size = index.size();
    std::vector<ScalarType*> values(size);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      values[v] = data_[index[v]];
    
    return(new MyMultiVec(Length_, values));
  }
  
  int GetVecLength () const
  {
    return(Length_);
  }
  
  int GetNumberVecs () const
  {
    return(NumberVecs_);
  }
  
  // Update *this with alpha * A * B + beta * (*this). 
  void MvTimesMatAddMv (const ScalarType alpha, const Anasazi::MultiVec<ScalarType> &A, 
			const Teuchos::SerialDenseMatrix<int, ScalarType> &B, 
			const ScalarType beta)
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert(MyA!=NULL);
    
    assert (Length_ == A.GetVecLength());
    assert (B.numRows() == A.GetNumberVecs());
    assert (B.numCols() == NumberVecs_);
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      {
        for (int i = 0 ; i < Length_ ; ++i)
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
  
  // Replace *this with alpha * A + beta * B. 
  void MvAddMv (const ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
		const ScalarType beta,  const Anasazi::MultiVec<ScalarType>& B)
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    MyMultiVec* MyB;
    MyB = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(B)); 
    assert (MyB != 0);
    
    assert (NumberVecs_ == A.GetNumberVecs());
    assert (NumberVecs_ == B.GetNumberVecs());
    
    assert (Length_ == A.GetVecLength());
    assert (Length_ == B.GetVecLength());
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      for (int i = 0 ; i < Length_ ; ++i)
	(*this)(i, v) = alpha * (*MyA)(i, v) + beta * (*MyB)(i, v);
  }
  
  // Compute a dense matrix B through the matrix-matrix multiply alpha * A^T * (*this). 
  void MvTransMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
		  Teuchos::SerialDenseMatrix< int, ScalarType >& B) const
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetVecLength() == Length_);
    assert (NumberVecs_ == B.numCols());
    assert (A.GetNumberVecs() == B.numRows());
    
    for (int v = 0 ; v < A.GetNumberVecs() ; ++v)
      {
        for (int w = 0 ; w < NumberVecs_ ; ++w)
	  {
	    ScalarType value = 0.0;
	    for (int i = 0 ; i < Length_ ; ++i)
	      {
		value += (*MyA)(i, v) * (*this)(i, w);
	      }
	    B(v, w) = alpha * value;
	  }
      }
  }
  
  // Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^T*this[i] where A[i] is the i-th column of A. 
  void MvDot (const Anasazi::MultiVec<ScalarType>& A, std::vector<ScalarType>* b) const
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (NumberVecs_ == (int)b->size());
    assert (NumberVecs_ == A.GetNumberVecs());
    assert (Length_ == A.GetVecLength());
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      {
        ScalarType value = 0.0;
        for (int i = 0 ; i < Length_ ; ++i)
          value += (*this)(i, v) * (*MyA)(i, v);
        (*b)[v] = value;
      }
  }
  
  void MvNorm (std::vector<magnitudeType> *normvec) const
  {
    assert (normvec != 0);
    assert (NumberVecs_ == (int)normvec->size());
    
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      {
        magnitudeType value = Teuchos::ScalarTraits<magnitudeType>::zero();
        for (int i = 0 ; i < Length_ ; ++i)
	  {
	    magnitudeType val = Teuchos::ScalarTraits<ScalarType>::magnitude((*this)(i, v));
	    value += val * val;
	  }
        (*normvec)[v] = Teuchos::ScalarTraits<magnitudeType>::squareroot(value);
      }
  }
  
  // Copy the vectors in A to a set of vectors in *this. The numvecs vectors in A are copied to a subset of vectors in *this indicated by the indices given in index.
  // FIXME: not so clear what the size of A and index.size() are...
  void SetBlock (const Anasazi::MultiVec<ScalarType>& A, 
		 const std::vector<int> &index)
  {
    MyMultiVec* MyA;
    MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetNumberVecs() >= (int)index.size());
    assert (A.GetVecLength() == Length_);
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
      for (int i = 0 ; i < Length_ ; ++i)
	(*this)(i, index[v])  = (*MyA)(i, v);
  }
  
  // Fill the vectors in *this with random numbers.
  virtual void  MvRandom ()
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      for (int i = 0 ; i < Length_ ; ++i)
	(*this)(i, v) = Teuchos::ScalarTraits<ScalarType>::random();
  }
  
  // Replace each element of the vectors in *this with alpha.
  virtual void  MvInit (const ScalarType alpha)
  {
    for (int v = 0 ; v < NumberVecs_ ; ++v)
      for (int i = 0 ; i < Length_ ; ++i)
	(*this)(i, v) = alpha;
  }
  
  void MvPrint (ostream &os) const
  {
    cout << "Number of rows = " << Length_ << endl;
    cout << "Number of vecs = " << NumberVecs_ << endl;
    
    for (int i = 0 ; i < Length_ ; ++i)
      {
        for (int v = 0 ; v < NumberVecs_ ; ++v)
          cout << (*this)(i, v) << " ";
        cout << endl;
      }
  }
  
  inline ScalarType& operator()(const int i, const int j)
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    
    return(data_[j][i]);
  }
  
  inline const ScalarType& operator()(const int i, const int j) const
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    
    return(data_[j][i]);
  }
  
  ScalarType* operator[](int v)
  {
    return(data_[v]);
  }
  
  ScalarType* operator[](int v) const
  {
    return(data_[v]);
  }
  
private:
  const int Length_;
  const int NumberVecs_;
  std::vector<ScalarType*> data_;
  std::vector<bool> ownership_;
};


#endif // MY_MULTIVECTOR_HPP
