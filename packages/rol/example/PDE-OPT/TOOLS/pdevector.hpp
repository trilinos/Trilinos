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

#ifndef ROL_PDEOPT_PDEVECTOR_HPP
#define ROL_PDEOPT_PDEVECTOR_HPP

#include "ROL_TpetraMultiVector.hpp"

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalSimVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualSimVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_DualSimVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_PrimalSimVector() {}

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap1(*pde);
    }

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDERieszMap1(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDERieszMap1(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class PDE_DualSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_PrimalSimVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_DualSimVector() {}

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap1(*pde);
    }

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDEInverseRieszMap1(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDEInverseRieszMap1(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

}; // class DualScaledTpetraMultiVector

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalOptVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualOptVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_DualOptVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_PrimalOptVector() {}

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap2(*pde);
    }

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDERieszMap2(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDERieszMap2(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class PDE_DualOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_PrimalOptVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_DualOptVector() {}

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap2(*pde);
    }

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDEInverseRieszMap2(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDEInverseRieszMap2(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

}; // class DualScaledTpetraMultiVector

#endif
