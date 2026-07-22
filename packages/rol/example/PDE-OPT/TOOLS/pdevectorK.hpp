// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_PDEVECTORK_HPP
#define ROL_PDEOPT_PDEVECTORK_HPP

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_PartitionedVector.hpp"

#include "assemblerK.hpp"
#include "solver.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double, Kokkos::HostSpace>;

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> PDEVectorRieszConstruct = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Riesz Construction Time");
    ROL::Ptr<Teuchos::Time> PDEVectorRieszApply     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Riesz Application Time");
    ROL::Ptr<Teuchos::Time> PDEVectorRieszSolve     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Riesz Solver Solution Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Construction Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Application Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorSimRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Solver Solution Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Construction Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Application Time");
//    ROL::Ptr<Teuchos::Time> PDEVectorOptRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Solver Solution Time");
  }
}
#endif

template <class Real, class DeviceType,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalVector;

template <class Real, class DeviceType,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualVector;

template <class Real, class DeviceType, class LO, class GO, class Node>
class PDE_PrimalVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    ROL::Ptr<Tpetra::CrsMatrix<>> RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<>> lumpedRiesz_;
    ROL::Ptr<Solver<Real>> solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<Tpetra::MultiVector<>> wex_;
    mutable bool isRieszVecInit_;

    mutable ROL::Ptr<PDE_DualVector<Real,DeviceType>> dual_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),1);
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    const Tpetra::MultiVector<>& applyRiesz(const Tpetra::MultiVector<> &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszApply);
      #endif
      if ( useRiesz_ ) {
	if (!isRieszVecInit_) {
          size_t n = in.getNumVectors();
          wex_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(
                 ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(), n);
	  isRieszVecInit_ = true;
	}
        if (useLumpedRiesz_) {
          wex_->elementWiseMultiply(static_cast<Real>(1), *(lumpedRiesz_->getVector(0)), in, static_cast<Real>(0));
        }
        else {
          RieszMap_->apply(in,*wex_);
        }
	return *wex_;
      }
      else {
        return in;
      }
    }

  public:
    virtual ~PDE_PrimalVector() {}

    using Assembler_type = Assembler<Real,DeviceType>;
    using PDE_type = PDE<Real,DeviceType>;
    using DynamicPDE_type = DynamicPDE<Real,DeviceType>;

    PDE_PrimalVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                     const ROL::Ptr<PDE_type> &pde,
                     const ROL::Ptr<Assembler_type> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isRieszVecInit_(false), isDualInitialized_(false) {}

    PDE_PrimalVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                     const ROL::Ptr<PDE_type> &pde,
                     const ROL::Ptr<Assembler_type> &assembler,
                     ROL::ParameterList &parlist,
                     bool sim = true)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), isRieszVecInit_(false), isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Lump Riesz Map", false);
      if (useRiesz_) {
        if (sim) assembler->assemblePDERieszMap1(RieszMap_, pde);
        else     assembler->assemblePDERieszMap2(RieszMap_, pde);
      }
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

    PDE_PrimalVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                     const ROL::Ptr<DynamicPDE_type> &pde,
                     Assembler_type &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isRieszVecInit_(false), isDualInitialized_(false) {}

    PDE_PrimalVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                     const ROL::Ptr<DynamicPDE_type> &pde,
                     Assembler_type &assembler,
                     ROL::ParameterList &parlist,
                     bool sim = true)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Lump Riesz Map", false);
      if (useRiesz_) {
        if (sim) assembler->assembleDynPDERieszMap1(RieszMap_, pde);
        else     assembler->assembleDynPDERieszMap2(RieszMap_, pde);
      }
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

    PDE_PrimalVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                     const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                     const ROL::Ptr<Solver<Real>> &solver,
                     const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isRieszVecInit_(false), isDualInitialized_(false) {
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
      const auto &ex = *dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const auto &ey = *ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector();
      size_t n = ex.getNumVectors();
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ex.dot( applyRiesz(ey), val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (const auto& v : val) xy += v;
      return xy;
    }

    ROL::Ptr<ROL::Vector<Real>> clone() const {
      auto &ey = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return ROL::makePtr<PDE_PrimalVector<Real,DeviceType,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = ROL::makePtr<PDE_DualVector<Real,DeviceType,LO,GO,Node>>(
                    ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),n),
                    RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      dual_vec_->getVector()->scale(static_cast<Real>(1),
        applyRiesz(*ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()));
      return *dual_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const auto &ex = dynamic_cast<const PDE_DualVector<Real,DeviceType,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_PrimalVector

template <class Real, class DeviceType, class LO, class GO, class Node>
class PDE_DualVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    ROL::Ptr<Tpetra::CrsMatrix<Real>> RieszMap_;
    ROL::Ptr<Tpetra::MultiVector<>> lumpedRiesz_;
    ROL::Ptr<Tpetra::MultiVector<>> recipLumpedRiesz_;
    ROL::Ptr<Solver<Real>> solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable ROL::Ptr<Tpetra::MultiVector<>> wex_;
    mutable bool isRieszVecInit_;

    mutable ROL::Ptr<PDE_PrimalVector<Real,DeviceType>> primal_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = ROL::makePtr<Tpetra::Vector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap());
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void invertLumpedRiesz(void) {
      recipLumpedRiesz_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),1);
      recipLumpedRiesz_->reciprocal(*lumpedRiesz_);
    }

    const ROL::Ptr<const Tpetra::MultiVector<>> applyRiesz(const ROL::Ptr<const Tpetra::MultiVector<>> &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszSolve);
      #endif
      if ( useRiesz_ ) {
	if (!isRieszVecInit_) {
          size_t n = in->getNumVectors();
          wex_ = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(
                 ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(), n);
	  isRieszVecInit_ = true;
	}
        if (useLumpedRiesz_) {
          wex_->elementWiseMultiply(static_cast<Real>(1), *(recipLumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          solver_->solve(wex_,in,false);
        }
	return wex_;
      }
      else {
        return in;
      }
    }

  public:
    virtual ~PDE_DualVector() {}

    using Assembler_type = Assembler<Real,DeviceType>;
    using PDE_type = PDE<Real,DeviceType>;
    using DynamicPDE_type = DynamicPDE<Real,DeviceType>;

    PDE_DualVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                   const ROL::Ptr<PDE_type> &pde,
                   const ROL::Ptr<Assembler_type> &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isRieszVecInit_(false), isDualInitialized_(false) {}

    PDE_DualVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                   const ROL::Ptr<PDE_type> &pde,
                   const ROL::Ptr<Assembler_type> &assembler,
                   ROL::ParameterList &parlist,
                   bool sim = true)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Lump Riesz Map", false);
      if (useRiesz_) {
        if (sim) assembler->assemblePDERieszMap1(RieszMap_, pde);
        else     assembler->assemblePDERieszMap2(RieszMap_, pde);
      }
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

    PDE_DualVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                   const ROL::Ptr<DynamicPDE_type> &pde,
                   Assembler_type &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(ROL::nullPtr),
        useRiesz_(false), useLumpedRiesz_(false), isRieszVecInit_(false), isDualInitialized_(false) {}

    PDE_DualVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                   const ROL::Ptr<DynamicPDE_type> &pde,
                   Assembler_type &assembler,
                   ROL::ParameterList &parlist,
                   bool sim = true)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), isRieszVecInit_(false),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist(sim ? "Sim" : "Opt").get("Lump Riesz Map", false);
      if (useRiesz_) {
        if (sim) assembler->assembleDynPDERieszMap1(RieszMap_, pde);
        else     assembler->assembleDynPDERieszMap2(RieszMap_, pde);
      }
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

    PDE_DualVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                   const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                   const ROL::Ptr<Solver<Real>> &solver,
                   const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isRieszVecInit_(false), isDualInitialized_(false) {
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
      const auto &ex = *dynamic_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector();
      const auto &ey = (ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey->getNumVectors();
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ex.dot( *applyRiesz(ey), val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (const auto& v : val) xy += v;
      return xy;
    }

    ROL::Ptr<ROL::Vector<Real>> clone() const {
      const auto &ey = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return ROL::makePtr<PDE_DualVector<Real,DeviceType,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),n),
             RieszMap_, solver_, lumpedRiesz_);
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = ROL::makePtr<PDE_PrimalVector<Real,DeviceType,LO,GO,Node>>(
                      ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(ROL::TpetraMultiVector<Real,LO,GO,Node>::getMap(),n),
                      RieszMap_, solver_, lumpedRiesz_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      primal_vec_->getVector()->scale(static_cast<Real>(1),
        *applyRiesz(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()));
      return *primal_vec_;
    }

    Real apply(const ROL::Vector<Real> &x) const {
      const auto &ex = dynamic_cast<const PDE_PrimalVector<Real,DeviceType,LO,GO,Node>&>(x);
      return ROL::TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }
}; // class PDE_DualVector

//using PDE_PrimalSimVector = PDE_PrimalVector<Real,DeviceType,LO,GO,Node>;
template <class Real, class DeviceType,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalSimVector : public PDE_PrimalVector<Real,DeviceType,LO,GO,Node> {
public:
  PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                      const ROL::Ptr<Assembler<Real,DeviceType>> &assembler)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                      const ROL::Ptr<Assembler<Real,DeviceType>> &assembler,
                      ROL::ParameterList &parlist)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,true) {}

  PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                      Assembler<Real,DeviceType> &assembler)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                      Assembler<Real,DeviceType> &assembler,
                      ROL::ParameterList &parlist)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,true) {}

  PDE_PrimalSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                      const ROL::Ptr<Solver<Real>> &solver,
                      const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,RieszMap,solver,lumpedRiesz) {}
};

//using PDE_DualSimVector = PDE_DualVector<Real,DeviceType,LO,GO,Node>;
template<class Real, class DeviceType,
         class LO=Tpetra::Map<>::local_ordinal_type, 
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type>
class PDE_DualSimVector : public PDE_DualVector<Real,DeviceType,LO,GO,Node> {
public:
  PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                    const ROL::Ptr<Assembler<Real,DeviceType>> &assembler)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                    const ROL::Ptr<Assembler<Real,DeviceType>> &assembler,
                    ROL::ParameterList &parlist)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,true) {}

  PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                    Assembler<Real,DeviceType> &assembler)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                    Assembler<Real,DeviceType> &assembler,
                    ROL::ParameterList &parlist)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,true) {}

  PDE_DualSimVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                    const ROL::Ptr<Solver<Real>> &solver,
                    const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,RieszMap,solver,lumpedRiesz) {}
};

//using PDE_PrimalOptVector = PDE_PrimalVector<Real,DeviceType,LO,GO,Node>;
template<class Real, class DeviceType,
         class LO=Tpetra::Map<>::local_ordinal_type, 
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type>
class PDE_PrimalOptVector : public PDE_PrimalVector<Real,DeviceType,LO,GO,Node> {
public:
  PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                      const ROL::Ptr<Assembler<Real,DeviceType>> &assembler)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                      const ROL::Ptr<Assembler<Real,DeviceType>> &assembler,
                      ROL::ParameterList &parlist)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,false) {}

  PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                      Assembler<Real,DeviceType> &assembler)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                      Assembler<Real,DeviceType> &assembler,
                      ROL::ParameterList &parlist)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,false) {}

  PDE_PrimalOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                      const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                      const ROL::Ptr<Solver<Real>> &solver,
                      const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
    : PDE_PrimalVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,RieszMap,solver,lumpedRiesz) {}
};

//using PDE_DualOptVector = PDE_DualVector<Real,DeviceType,LO,GO,Node>;
template <class Real, class DeviceType,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualOptVector : public PDE_DualVector<Real,DeviceType,LO,GO,Node> {
public:
  PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                    const ROL::Ptr<Assembler<Real,DeviceType>> &assembler)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                    const ROL::Ptr<Assembler<Real,DeviceType>> &assembler,
                    ROL::ParameterList &parlist)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,false) {}

  PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                    Assembler<Real,DeviceType> &assembler)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler) {}

  PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<DynamicPDE<Real,DeviceType>> &pde,
                    Assembler<Real,DeviceType> &assembler,
                    ROL::ParameterList &parlist)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,pde,assembler,parlist,false) {}

  PDE_DualOptVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &tpetra_vec,
                    const ROL::Ptr<Tpetra::CrsMatrix<>> &RieszMap,
                    const ROL::Ptr<Solver<Real>> &solver,
                    const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node>> &lumpedRiesz)
    : PDE_DualVector<Real,DeviceType,LO,GO,Node>(tpetra_vec,RieszMap,solver,lumpedRiesz) {}
};

template<class Real,
         class LO=Tpetra::Map<>::local_ordinal_type, 
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type>
class PDE_OptVector : public ROL::PartitionedVector<Real> {
private:
  const int rank_;

  mutable ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node>> dvec1_;
  mutable ROL::Ptr<ROL::StdVector<Real>> dvec2_;
  mutable ROL::Ptr<PDE_OptVector<Real,LO,GO,Node>> dvec_;
  mutable bool isDualInit_;

public:
  PDE_OptVector(const ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node>> &vec1,
                const ROL::Ptr<ROL::StdVector<Real>>                    &vec2,
                int rank = 0 ) 
    : ROL::PartitionedVector<Real>({vec1,vec2}), rank_(rank), isDualInit_(false) {}

  ROL::Ptr<const ROL::TpetraMultiVector<Real,LO,GO,Node>> getField(void) const { 
    return ROL::staticPtrCast<const ROL::TpetraMultiVector<Real,LO,GO,Node>>(ROL::PartitionedVector<Real>::get(0));
  }

  ROL::Ptr<const ROL::StdVector<Real>> getParameter(void) const { 
    return ROL::staticPtrCast<const ROL::StdVector<Real>>(ROL::PartitionedVector<Real>::get(1));
  }

  ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node>> getField(void) { 
    return ROL::staticPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node>>(ROL::PartitionedVector<Real>::get(0));
  }

  ROL::Ptr<ROL::StdVector<Real>> getParameter(void) { 
    return ROL::staticPtrCast<ROL::StdVector<Real>>(ROL::PartitionedVector<Real>::get(1));
  }

  ROL::Ptr<ROL::Vector<Real>> clone(void) const {
    auto pvec1 = ROL::staticPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node>>(getField()->clone());
    auto pvec2 = ROL::staticPtrCast<ROL::StdVector<Real>>(getParameter()->clone());
    return ROL::makePtr<PDE_OptVector<Real,LO,GO,Node>>(pvec1,pvec2,rank_);
  }

  const ROL::Vector<Real>& dual(void) const {
    if (!isDualInit_) {
      dvec1_ = ROL::staticPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node>>(getField()->dual().clone());
      dvec2_ = ROL::staticPtrCast<ROL::StdVector<Real>>(getParameter()->dual().clone());
      dvec_  = ROL::makePtr<PDE_OptVector<Real,LO,GO,Node>>(dvec1_,dvec2_,rank_);
      isDualInit_ = true;
    }
    dvec1_->set(getField()->dual());
    dvec2_->set(getParameter()->dual());
    return *dvec_;
  }

  void setField(const ROL::Vector<Real>& vec) {
    auto vec_cast = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real,LO,GO,Node>>(vec);
    ROL::PartitionedVector<Real>::set(0,vec_cast);
  }
  
  void setParameter(const ROL::Vector<Real>& vec) { 
    auto vec_cast = ROL::dynamicPtrCast<ROL::StdVector<Real>>(vec);
    ROL::PartitionedVector<Real>::set(1,vec_cast);
  }

  void print(std::ostream &outStream) const {
    getField()->print(outStream);
    if (rank_==0) getParameter()->print(outStream);
  }
}; // class PDE_OptVector

#endif
