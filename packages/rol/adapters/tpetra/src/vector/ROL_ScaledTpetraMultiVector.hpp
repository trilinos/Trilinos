// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_SCALEDTPETRAMULTIVECTOR_HPP
#define ROL_SCALEDTPETRAMULTIVECTOR_HPP

/** \class ROL::ScaledTpetraMultiVector
    \brief Implements the ROL::Vector interface for a Tpetra_MultiVector
           with scaled dot function.
    \author Created by Drew Kouri
*/

#include "ROL_TpetraMultiVector.hpp"

namespace ROL {

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PrimalScaledTpetraMultiVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class DualScaledTpetraMultiVector;

template <class Real, class LO, class GO, class Node>
class PrimalScaledTpetraMultiVector : public TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<const Tpetra::Vector<Real,LO,GO,Node> > scale_vec_;
    mutable Teuchos::RCP<DualScaledTpetraMultiVector<Real> > dual_vec_;

  public:

    PrimalScaledTpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                                  const Teuchos::RCP<const Tpetra::Vector<Real,LO,GO,Node> > &scale_vec)
      : TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), scale_vec_(scale_vec) {}

    Real dot( const Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Tpetra::MultiVector<Real,LO,GO,Node> &ex
        = *(Teuchos::dyn_cast<const TpetraMultiVector<Real,LO,GO,Node> >(x).getVector());
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      // Scale x with scale_vec_
      Tpetra::MultiVector<Real,LO,GO,Node> wex(ex);
      wex.elementWiseMultiply(1,*scale_vec_,ex,0);
      // Perform Euclidean dot between *this and scaled x
      Teuchos::Array<Real> val(1, 0);
      ey.dot( wex, val );
      return val[0];
    }

    Teuchos::RCP<Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PrimalScaledTpetraMultiVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ey.getMap(),n,false)),
             scale_vec_));
    }

    const Vector<Real> & dual() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      // Scale *this with scale_vec_
      Tpetra::MultiVector<Real,LO,GO,Node> wey(ey);
      wey.elementWiseMultiply(1,*scale_vec_,ey,0);
      // Return a copy of scaled *this
      dual_vec_ = Teuchos::rcp(new DualScaledTpetraMultiVector<Real,LO,GO,Node>(
                  Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(wey)),
                  scale_vec_));
      return *dual_vec_;
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class DualScaledTpetraMultiVector : public TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<const Tpetra::Vector<Real,LO,GO,Node> > scale_vec_;
    mutable Teuchos::RCP<PrimalScaledTpetraMultiVector<Real> > primal_vec_;

  public:

    DualScaledTpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                                const Teuchos::RCP<const Tpetra::Vector<Real,LO,GO,Node> > &scale_vec)
      : TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), scale_vec_(scale_vec) {}

    Real dot( const Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Tpetra::MultiVector<Real,LO,GO,Node> &ex
        = *(Teuchos::dyn_cast<const TpetraMultiVector<Real,LO,GO,Node> >(x).getVector());
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      // Scale x with 1/scale_vec_
      Tpetra::Vector<Real,LO,GO,Node> W(*scale_vec_);
      W.reciprocal(*scale_vec_);
      Tpetra::MultiVector<Real,LO,GO,Node> wex(ex);
      wex.elementWiseMultiply(1,W,ex,0);
      // Perform Euclidean dot between *this and scaled x
      Teuchos::Array<Real> val(1, 0);
      ey.dot( wex, val );
      return val[0];
    }

    Teuchos::RCP<Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new DualScaledTpetraMultiVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ey.getMap(),n,false)),
             scale_vec_));
    }

    const Vector<Real> & dual() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      // Scale *this with 1/scale_vec_
      Tpetra::Vector<Real,LO,GO,Node> W(*scale_vec_);
      W.reciprocal(*scale_vec_);
      Tpetra::MultiVector<Real,LO,GO,Node> wey(ey);
      wey.elementWiseMultiply(1,W,ey,0);
      // Return a copy of scaled *this
      primal_vec_ = Teuchos::rcp(new PrimalScaledTpetraMultiVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(wey)),
                    scale_vec_));
      return *primal_vec_;
    }

}; // class DualScaledTpetraMultiVector

} // namespace ROL

#endif
