#ifndef TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
#define TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP

#include <Teuchos_ScalarTraits.hpp>

template<class Scalar, class LO>
class LittleBlock {
private:
  typedef Teuchos::ScalarTraits<Scalar> STS;

public:
  LittleBlock (Scalar* const A, const LO blockSize, const LO strideX, const LO strideY) :
    A_ (A), blockSize_ (blockSize), strideX_ (strideX), strideY_ (strideY) {}

  Scalar* getRawPtr () const {
    return A_;
  }

  Scalar& operator() (const LO i, const LO j) const {
    return A_[i * strideX_ + j * strideY_];
  }

  template<class LittleBlockType>
  void update (const Scalar& alpha, const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) += alpha * X(i,j);
      }
    }
  }

  template<class LittleBlockType>
  void assign (const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) = X(i,j);
      }
    }
  }

  void scale (const Scalar& alpha) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) *= alpha;
      }
    }
  }

  void fill (const Scalar& alpha) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) = alpha;
      }
    }
  }

  template<class LittleBlockType>
  void absmax (const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        Scalar& Y_ij = (*this)(i,j);
        const Scalar X_ij = X(i,j);
        Y_ij = std::max (STS::magnitude (Y_ij), STS::magnitude (X_ij));
      }
    }
  }

private:
  Scalar* const A_;
  const LO blockSize_;
  const LO strideX_;
  const LO strideY_;
};



template<class Scalar, class LO>
class LittleVector {
private:
  typedef Teuchos::ScalarTraits<Scalar> STS;

public:
  LittleVector (Scalar* const A, const LO blockSize, const LO strideX) :
    A_ (A), blockSize_ (blockSize), strideX_ (strideX) {}

  Scalar* getRawPtr () const {
    return A_;
  }

  Scalar& operator() (const LO i) const {
    return A_[i * strideX_];
  }

  template<class LittleBlockType>
  void update (const Scalar& alpha, const LittleBlockType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) += alpha * X(i);
    }
  }

  template<class LittleBlockType>
  void assign (const LittleBlockType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = X(i);
    }
  }

  void scale (const Scalar& alpha) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) *= alpha;
    }
  }

  void fill (const Scalar& alpha) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = alpha;
    }
  }

  template<class LittleBlockType>
  void absmax (const LittleBlockType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      Scalar& Y_i = (*this)(i);
      Y_i = std::max (STS::magnitude (Y_i), STS::magnitude (X (i)));
    }
  }

private:
  Scalar* const A_;
  const LO blockSize_;
  const LO strideX_;
};


#endif // TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
