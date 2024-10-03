// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_PDEVECTOR_HPP
#define ROL_PDEOPT_PDEVECTOR_HPP

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"

#include "assembler.hpp"
#include "solver.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Construction Time");
    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Application Time");
    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Solver Solution Time");
    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Construction Time");
    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Application Time");
    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Solver Solution Time");
  }
}
#endif


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
    ROL::Ptr<Tpetra::CrsMatrix<> > RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<> > lumpedRiesz_;
    ROL::Ptr<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<PDE_DualSimVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),1);
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void applyRiesz(const ROL::Ptr<Tpetra::MultiVector<> > &out,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszApply);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(lumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          RieszMap_->apply(*in,*out);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_PrimalSimVector() {}

    PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<PDE<Real> > &pde,
                        const ROL::Ptr<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<PDE<Real> > &pde,
                        const ROL::Ptr<Assembler<Real> > &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      if (useRiesz_) assembler->assemblePDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<DynamicPDE<Real> > &pde,
                        Assembler<Real> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<DynamicPDE<Real> > &pde,
                        Assembler<Real> &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      if (useRiesz_) assembler.assembleDynPDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<Tpetra::CrsMatrix<> > &RieszMap,
                        const ROL::Ptr<Solver<Real> > &solver,
                        const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != ROL::nullPtr) {
        useLumpedRiesz_ = (lumpedRiesz_ != ROL::nullPtr);
        useRiesz_ = (solver_ != ROL::nullPtr) || useLumpedRiesz_;
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(), n);
      applyRiesz(wex,ex);
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

    ROL::Ptr<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return ROL::makePtr<PDE_PrimalSimVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = ROL::makePtr<PDE_DualSimVector<Real,LO,GO,Node>>(
                    ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
                    RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const PDE_DualSimVector<Real,LO,GO,Node> &ex = dynamic_cast<const PDE_DualSimVector<Real,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_PrimalSimVector

template <class Real, class LO, class GO, class Node>
class PDE_DualSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    ROL::Ptr<Tpetra::CrsMatrix<Real> > RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<> > lumpedRiesz_;
    ROL::Ptr<Tpetra::MultiVector<> > recipLumpedRiesz_;
    ROL::Ptr<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<PDE_PrimalSimVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::Vector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap());
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void invertLumpedRiesz(void) {
      recipLumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),1);
      recipLumpedRiesz_->reciprocal(*lumpedRiesz_);
    }

    void applyRiesz(const ROL::Ptr<Tpetra::MultiVector<> > &out,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszSolve);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(recipLumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          solver_->solve(out,in,false);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_DualSimVector() {}

    PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<PDE<Real> > &pde,
                      const ROL::Ptr<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<PDE<Real> > &pde,
                      const ROL::Ptr<Assembler<Real> > &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      if (useRiesz_) assembler->assemblePDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real> > &pde,
                      Assembler<Real> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real> > &pde,
                      Assembler<Real> &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      assembler.assembleDynPDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<Tpetra::CrsMatrix<> > &RieszMap,
                      const ROL::Ptr<Solver<Real> > &solver,
                      const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != ROL::nullPtr) {
        useLumpedRiesz_ = (lumpedRiesz_ != ROL::nullPtr);
        useRiesz_ = (solver_ != ROL::nullPtr) || useLumpedRiesz_;
        if (useLumpedRiesz_) {
          invertLumpedRiesz();
        }
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(), n);
      applyRiesz(wex,ex);
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

    ROL::Ptr<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return ROL::makePtr<PDE_DualSimVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = ROL::makePtr<PDE_PrimalSimVector<Real,LO,GO,Node>>(
                      ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
                      RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const PDE_PrimalSimVector<Real,LO,GO,Node> &ex = dynamic_cast<const PDE_PrimalSimVector<Real,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_DualSimVector

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
    ROL::Ptr<Tpetra::CrsMatrix<> > RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<> > lumpedRiesz_;
    ROL::Ptr<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<PDE_DualOptVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),1);
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void applyRiesz(const ROL::Ptr<Tpetra::MultiVector<> > &out,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszApply);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(lumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          RieszMap_->apply(*in,*out);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_PrimalOptVector() {}

    PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<PDE<Real> > &pde,
                        const ROL::Ptr<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<PDE<Real> > &pde,
                        const ROL::Ptr<Assembler<Real> > &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      if (useRiesz_) assembler->assemblePDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<DynamicPDE<Real> > &pde,
                        Assembler<Real> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<DynamicPDE<Real> > &pde,
                        Assembler<Real> &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      assembler.assembleDynPDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const ROL::Ptr<Tpetra::CrsMatrix<> > &RieszMap,
                        const ROL::Ptr<Solver<Real> > &solver,
                        const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != ROL::nullPtr) {
        useLumpedRiesz_ = (lumpedRiesz_ != ROL::nullPtr);
        useRiesz_ = (solver_ != ROL::nullPtr) || useLumpedRiesz_;
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(), n);
      applyRiesz(wex,ex);
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

    ROL::Ptr<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return ROL::makePtr<PDE_PrimalOptVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = ROL::makePtr<PDE_DualOptVector<Real,LO,GO,Node>>(
                    ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
                    RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const PDE_DualOptVector<Real,LO,GO,Node> &ex = dynamic_cast<const PDE_DualOptVector<Real,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_PrimalOptVector

template <class Real, class LO, class GO, class Node>
class PDE_DualOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    ROL::Ptr<Tpetra::CrsMatrix<Real> > RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<> > lumpedRiesz_;
    ROL::Ptr<Tpetra::MultiVector<> > recipLumpedRiesz_;
    ROL::Ptr<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<PDE_PrimalOptVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::Vector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap());
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void invertLumpedRiesz(void) {
      recipLumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),1);
      recipLumpedRiesz_->reciprocal(*lumpedRiesz_);
    }

    void applyRiesz(const ROL::Ptr<Tpetra::MultiVector<> > &out,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszSolve);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(recipLumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          solver_->solve(out,in,false);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_DualOptVector() {}

    PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<PDE<Real> > &pde,
                      const ROL::Ptr<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<PDE<Real> > &pde,
                      const ROL::Ptr<Assembler<Real> > &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      if (useRiesz_) assembler->assemblePDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real> > &pde,
                      Assembler<Real> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real> > &pde,
                      Assembler<Real> &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      assembler.assembleDynPDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != ROL::nullPtr);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const ROL::Ptr<Tpetra::CrsMatrix<> > &RieszMap,
                      const ROL::Ptr<Solver<Real> > &solver,
                      const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != ROL::nullPtr) {
        useLumpedRiesz_ = (lumpedRiesz_ != ROL::nullPtr);
        useRiesz_ = (solver_ != ROL::nullPtr) || useLumpedRiesz_;
        if (useLumpedRiesz_) {
          invertLumpedRiesz();
        }
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(), n);
      applyRiesz(wex,ex);
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

    ROL::Ptr<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return ROL::makePtr<PDE_DualOptVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = ROL::makePtr<PDE_PrimalOptVector<Real,LO,GO,Node>>(
                      ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real>::getMap(),n),
                      RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const PDE_PrimalOptVector<Real,LO,GO,Node> &ex = dynamic_cast<const PDE_PrimalOptVector<Real,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_DualOptVector

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_OptVector : public ROL::Vector<Real> {
private:
  ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > vec1_;
  ROL::Ptr<ROL::StdVector<Real> >                    vec2_;
  const int rank_;
  mutable ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > dual_vec1_;
  mutable ROL::Ptr<ROL::StdVector<Real> >                    dual_vec2_;
  mutable ROL::Ptr<PDE_OptVector<Real,LO,GO,Node> >          dual_vec_;
  mutable bool isDualInitialized_;

public:
  PDE_OptVector(const ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > &vec1,
                const ROL::Ptr<ROL::StdVector<Real> >                    &vec2,
                const int rank = 0 ) 
    : vec1_(vec1), vec2_(vec2), rank_(rank), isDualInitialized_(false) {

    dual_vec1_ = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->dual().clone());
    dual_vec2_ = ROL::dynamicPtrCast<ROL::StdVector<Real> >(vec2_->dual().clone());
  }

  PDE_OptVector(const ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > &vec)
    : vec1_(vec), vec2_(ROL::nullPtr), rank_(0), dual_vec2_(ROL::nullPtr), isDualInitialized_(false) {
    dual_vec1_ = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->dual().clone());
  }

  PDE_OptVector(const ROL::Ptr<ROL::StdVector<Real> > &vec,
                const int rank = 0)
    : vec1_(ROL::nullPtr), vec2_(vec), rank_(rank), dual_vec1_(ROL::nullPtr), isDualInitialized_(false) {
    dual_vec2_ = ROL::dynamicPtrCast<ROL::StdVector<Real> >(vec2_->dual().clone());
  }

  void set( const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->set(*(xs.getField()));
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->set(*(xs.getParameter()));
    }
  }

  void plus( const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->plus(*(xs.getField()));
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->plus(*(xs.getParameter()));
    }
  }

  void scale( const Real alpha ) {
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->scale(alpha);
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->scale(alpha);
    }
  }

  void axpy( const Real alpha, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->axpy(alpha,*(xs.getField()));
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->axpy(alpha,*(xs.getParameter()));
    }
  }

  Real dot( const ROL::Vector<Real> &x ) const {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    Real val(0);
    if ( vec1_ != ROL::nullPtr ) {
      val += vec1_->dot(*(xs.getField()));
    }
    if ( vec2_ != ROL::nullPtr ) {
      val += vec2_->dot(*(xs.getParameter()));
    }
    return val;
  }

  Real norm() const {
    Real val(0);
    if ( vec1_ != ROL::nullPtr ) {
      Real norm1 = vec1_->norm();
      val += norm1*norm1;
    }
    if ( vec2_ != ROL::nullPtr ) {
      Real norm2 = vec2_->norm();
      val += norm2*norm2;
    }
    return std::sqrt(val);
  } 

  ROL::Ptr<ROL::Vector<Real> > clone(void) const {
    if ( vec2_ == ROL::nullPtr ) {
      return ROL::makePtr<PDE_OptVector<Real,LO,GO,Node>>(
             ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->clone()));
    }
    if ( vec1_ == ROL::nullPtr ) {
      return ROL::makePtr<PDE_OptVector<Real,LO,GO,Node>>(
             ROL::dynamicPtrCast<ROL::StdVector<Real> >(vec2_->clone()));
    }
    return ROL::makePtr<PDE_OptVector<Real,LO,GO,Node>>(
           ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->clone()),
           ROL::dynamicPtrCast<ROL::StdVector<Real> >(vec2_->clone()));
  }

  const ROL::Vector<Real> & dual(void) const {
    if ( !isDualInitialized_ ) {
      if ( vec1_ == ROL::nullPtr ) {
        dual_vec_ = ROL::makePtr<PDE_OptVector<Real>>(dual_vec2_);
      }
      else if ( vec2_ == ROL::nullPtr ) {
        dual_vec_ = ROL::makePtr<PDE_OptVector<Real>>(dual_vec1_);
      }
      else {
        dual_vec_ = ROL::makePtr<PDE_OptVector<Real>>(dual_vec1_,dual_vec2_);
      }
      isDualInitialized_ = true;
    }
    if ( vec1_ != ROL::nullPtr ) {
      dual_vec1_->set(vec1_->dual());
    }
    if ( vec2_ != ROL::nullPtr ) {
      dual_vec2_->set(vec2_->dual());
    }
    return *dual_vec_;
  }

  Real apply(const ROL::Vector<Real> &x) const {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    Real val(0);
    if ( vec1_ != ROL::nullPtr ) val += vec1_->apply(*(xs.getField()));
    if ( vec2_ != ROL::nullPtr ) val += vec2_->apply(*(xs.getParameter()));
    return val;
  }

  ROL::Ptr<ROL::Vector<Real> > basis( const int i )  const {
    ROL::Ptr<ROL::Vector<Real> > e;
    if ( vec1_ != ROL::nullPtr && vec2_ != ROL::nullPtr ) {
      int n1 = vec1_->dimension();
      ROL::Ptr<ROL::Vector<Real> > e1, e2;
      if ( i < n1 ) {
        e1 = vec1_->basis(i);
        e2 = vec2_->clone(); e2->zero();
      }
      else {
        e1 = vec1_->clone(); e1->zero();
        e2 = vec2_->basis(i-n1);
      }
      e = ROL::makePtr<PDE_OptVector>(
        ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real> >(e1),
        ROL::dynamicPtrCast<ROL::StdVector<Real> >(e2));
    }
    if ( vec1_ != ROL::nullPtr && vec2_ == ROL::nullPtr ) {
      int n1 = vec1_->dimension();
      ROL::Ptr<ROL::Vector<Real> > e1;
      if ( i < n1 ) {
        e1 = vec1_->basis(i);
      }
      else {
        e1->zero();
      }
      e = ROL::makePtr<PDE_OptVector>(
        ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real> >(e1));
    }
    if ( vec1_ == ROL::nullPtr && vec2_ != ROL::nullPtr ) {
      int n2 = vec2_->dimension();
      ROL::Ptr<ROL::Vector<Real> > e2;
      if ( i < n2 ) {
        e2 = vec2_->basis(i);
      }
      else {
        e2->zero();
      }
      e = ROL::makePtr<PDE_OptVector>(
        ROL::dynamicPtrCast<ROL::StdVector<Real> >(e2));
    }
    return e;
  }

  void applyUnary( const ROL::Elementwise::UnaryFunction<Real> &f ) {
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->applyUnary(f);
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->applyUnary(f);
    }
  }

  void applyBinary( const ROL::Elementwise::BinaryFunction<Real> &f, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = dynamic_cast<const PDE_OptVector<Real>&>(x);
    if ( vec1_ != ROL::nullPtr ) {
      vec1_->applyBinary(f,*xs.getField());
    }
    if ( vec2_ != ROL::nullPtr ) {
      vec2_->applyBinary(f,*xs.getParameter());
    }
  }

  Real reduce( const ROL::Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    if ( vec1_ != ROL::nullPtr ) {
      r.reduce(vec1_->reduce(r),result);
    }
    if ( vec2_ != ROL::nullPtr ) {
      r.reduce(vec2_->reduce(r),result);
    }
    return result;
  }

  int dimension() const {
    int dim(0);
    if ( vec1_ != ROL::nullPtr ) {
      dim += vec1_->dimension();
    }
    if ( vec2_ != ROL::nullPtr ) {
      dim += vec2_->dimension();
    }
    return dim;
  }

  void randomize(const Real l = 0.0, const Real u = 1.0) {
    if (vec1_ != ROL::nullPtr) {
      vec1_->randomize(l,u);
    }
    if (vec2_ != ROL::nullPtr) {
      vec2_->randomize(l,u);
    }
  }

  void print(std::ostream &outStream) const {
    if (vec1_ != ROL::nullPtr) {
      vec1_->print(outStream);
    }
    if (vec2_ != ROL::nullPtr) {
      if (rank_ == 0) {
        vec2_->print(outStream);
      }
    }
  }

  ROL::Ptr<const ROL::TpetraMultiVector<Real> > getField(void) const { 
    return vec1_;
  }

  ROL::Ptr<const ROL::StdVector<Real> > getParameter(void) const { 
    return vec2_; 
  }

  ROL::Ptr<ROL::TpetraMultiVector<Real> > getField(void) { 
    return vec1_;
  }

  ROL::Ptr<ROL::StdVector<Real> > getParameter(void) { 
    return vec2_; 
  }

  void setField(const ROL::Vector<Real>& vec) { 
    vec1_->set(vec);
  }
  
  void setParameter(const ROL::Vector<Real>& vec) { 
    vec2_->set(vec); 
  }
}; // class PDE_OptVector

#endif
