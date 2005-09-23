#ifndef MYMULTIVECTOR_HPP
#define MYMULTIVECTOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBlockDavidson.hpp"

#include "AnasaziMultiVec.hpp"
//#include "AnasaziOperator.hpp"
//#include "AnasaziOperatorTraits.hpp"
//#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include <complex>

// Wrapper for Teuchos_SerialDenseMatrix<type>
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

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
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
      data_.resize(GetNumberVecs());
      ownership_.resize(GetNumberVecs());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        data_[v] = rhs[v];
        ownership_[v] = false;
      }
    }

    MyMultiVec(const MyMultiVec& rhs) :
      Length_(rhs.GetVecLength()),
      NumberVecs_(rhs.GetNumberVecs())
    {
      data_.resize(GetNumberVecs());
      ownership_.resize(GetNumberVecs());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        data_[v] = new ScalarType[GetVecLength()];
        ownership_[v] = true;
      }

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*this)(i, v) = rhs(i, v);
      }
    }

    ~MyMultiVec()
    {
      // FIXME!
//        for (int v = 0 ; v < GetNumberVecs() ; ++v)
  //        if (ownership_[v]) delete[] data_[v];
    }

    MyMultiVec* Clone(const int NumberVecs) const
    {
      MyMultiVec* tmp = new MyMultiVec(GetVecLength(), NumberVecs);

      //   for (int v = 0 ; v < NumberVecs ; ++v)
      //         for (int i = 0 ; i < GetVecLength() ; ++i)
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
      MyMultiVec* tmp = new MyMultiVec(GetVecLength(), size);

      for (int v = 0 ; v < index.size() ; ++v)
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*tmp)(i, v) = (*this)(i, index[v]);

      return(tmp);
    }

    MyMultiVec* CloneView(const std::vector< int > &index) 
    {
      int size = index.size();
      std::vector<ScalarType*> values(size);

      for (int v = 0 ; v < index.size() ; ++v)
        values[v] = data_[index[v]];

      return(new MyMultiVec(GetVecLength(), values));
    }

  const MyMultiVec* CloneView(const std::vector< int > &index) const
    {
      int size = index.size();
      std::vector<ScalarType*> values(size);

      for (int v = 0 ; v < index.size() ; ++v)
        values[v] = data_[index[v]];

      return(new MyMultiVec(GetVecLength(), values));
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

      assert (GetVecLength() == A.GetVecLength());
      assert (B.numRows() == A.GetNumberVecs());
      assert (B.numCols() == GetNumberVecs());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < GetVecLength() ; ++i)
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

      assert (GetNumberVecs() == A.GetNumberVecs());
      assert (GetNumberVecs() == B.GetNumberVecs());

      assert (GetVecLength() == A.GetVecLength());
      assert (GetVecLength() == B.GetVecLength());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*this)(i, v) = alpha * (*MyA)(i, v) + beta * (*MyB)(i, v);
    }

    // Compute a dense matrix B through the matrix-matrix multiply alpha * A^T * (*this). 
    void MvTransMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
                    Teuchos::SerialDenseMatrix< int, ScalarType >& B) const
    {
      MyMultiVec* MyA;
      MyA = dynamic_cast<MyMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
      assert (MyA != 0);

      assert (A.GetVecLength() == GetVecLength());
      assert (GetNumberVecs() == B.numCols());
      assert (A.GetNumberVecs() == B.numRows());

      for (int v = 0 ; v < A.GetNumberVecs() ; ++v)
      {
        for (int w = 0 ; w < GetNumberVecs() ; ++w)
        {
          ScalarType value = 0.0;
          for (int i = 0 ; i < GetVecLength() ; ++i)
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

      assert (GetNumberVecs() == b->size());
      assert (GetNumberVecs() == A.GetNumberVecs());
      assert (GetVecLength() == A.GetVecLength());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        ScalarType value = 0.0;
        for (int i = 0 ; i < GetVecLength() ; ++i)
          value += (*this)(i, v) * (*MyA)(i, v);
        (*b)[v] = value;
      }
    }

  void MvNorm (std::vector<magnitudeType> *normvec) const
    {
      assert (GetNumberVecs() == normvec->size());

      for (int v = 0 ; v < GetNumberVecs() ; ++v)
      {
        magnitudeType value = Teuchos::ScalarTraits<magnitudeType>::zero();
        for (int i = 0 ; i < GetVecLength() ; ++i)
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

      assert (A.GetNumberVecs() >= index.size());
      assert (A.GetVecLength() == GetVecLength());

      for (int v = 0 ; v < index.size() ; ++v)
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*this)(i, index[v])  = (*MyA)(i, v);
    }

    // Fill the vectors in *this with random numbers.
    virtual void  MvRandom ()
    {
      for (int v = 0 ; v < GetNumberVecs() ; ++v)
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*this)(i, v) = Teuchos::ScalarTraits<ScalarType>::random();
    }

    // Replace each element of the vectors in *this with alpha.
    virtual void  MvInit (const ScalarType alpha)
    {
      for (int v = 0 ; v < GetNumberVecs() ; ++v)
        for (int i = 0 ; i < GetVecLength() ; ++i)
          (*this)(i, v) = alpha;
    }

    void       MvPrint (ostream &os) const
    {
      cout << "Number of rows = " << GetVecLength() << endl;
      cout << "Number of vecs = " << GetNumberVecs() << endl;

      for (int i = 0 ; i < GetVecLength() ; ++i)
      {
        for (int v = 0 ; v < GetNumberVecs() ; ++v)
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






// ============================================================================ 

namespace Anasazi {

template<>
class MultiVecTraits<ScalarType, MyMultiVec>
{
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType magnitudeType;

  public:

    static Teuchos::RefCountPtr<MyMultiVec> Clone( const MyMultiVec& mv, const int numvecs )
    { return Teuchos::rcp( new MyMultiVec(mv.GetVecLength(), numvecs) ); }

    static Teuchos::RefCountPtr<MyMultiVec> CloneCopy( const MyMultiVec& mv )
    { return Teuchos::rcp( new MyMultiVec( mv ) ); }

    static Teuchos::RefCountPtr<MyMultiVec> CloneCopy( const MyMultiVec& mv, const std::vector<int>& index )
    { 
      return(Teuchos::rcp(mv.CloneCopy(index)));
    }

    static Teuchos::RefCountPtr<MyMultiVec> CloneView( MyMultiVec& mv, const std::vector<int>& index )
    { 
      return(Teuchos::rcp(mv.CloneView(index)));
    }

    static Teuchos::RefCountPtr<const MyMultiVec> CloneView( const MyMultiVec& mv, const std::vector<int>& index )
    { 
      return(Teuchos::rcp(mv.CloneView(index)));
    }

    static int GetVecLength( const MyMultiVec& mv )
    { 
      return(mv.GetVecLength());
    }

    static int GetNumberVecs( const MyMultiVec& mv )
    { 
      return(mv.GetNumberVecs());
    }

    static void MvTimesMatAddMv( const ScalarType alpha, const MyMultiVec& A, 
                                const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                                const ScalarType beta, MyMultiVec& mv )
    { 
      mv.MvTimesMatAddMv(alpha, A, B, beta);
    }

    static void MvAddMv( const ScalarType alpha, const MyMultiVec& A, const ScalarType beta, const MyMultiVec& B, MyMultiVec& mv )
    { 
      mv.MvAddMv(alpha, A, beta, B);
    }

    static void MvTransMv( const ScalarType alpha, const MyMultiVec& A, const MyMultiVec& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B )
    { 
      mv.MvTransMv(alpha, A, B);
    }

    static void MvDot( const MyMultiVec& mv, const MyMultiVec& A, std::vector<ScalarType>* b )
    {
      mv.MvDot(A, b);
    }
  
  static void MvNorm( const MyMultiVec& mv, std::vector<magnitudeType>* normvec )
  { 
      mv.MvNorm(normvec);
    }

    static void SetBlock( const MyMultiVec& A, const std::vector<int>& index, MyMultiVec& mv )
    { 
      mv.SetBlock(A, index);
    }

    static void MvRandom( MyMultiVec& mv )
    { 
      mv.MvRandom(); 
    }

    static void MvInit( MyMultiVec& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { 
      mv.MvInit(alpha);
    }

    static void MvPrint( const MyMultiVec& mv, ostream& os )
    { 
      mv.MvPrint(os); 
    }

};        

}


#endif // MYMULTIVECTOR_HPP
