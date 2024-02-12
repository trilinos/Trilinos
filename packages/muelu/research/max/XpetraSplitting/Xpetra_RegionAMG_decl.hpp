// @HEADER
// //
// // ***********************************************************************
// //
// //        MueLu: A package for multigrid based preconditioning
// //                  Copyright 2012 Sandia Corporation
// //
// // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// // the U.S. Government retains certain rights in this software.
// //
// // Redistribution and use in source and binary forms, with or without
// // modification, are permitted provided that the following conditions are
// // met:
// //
// // 1. Redistributions of source code must retain the above copyright
// // notice, this list of conditions and the following disclaimer.
// //
// // 2. Redistributions in binary form must reproduce the above copyright
// // notice, this list of conditions and the following disclaimer in the
// // documentation and/or other materials provided with the distribution.
// //
// // 3. Neither the name of the Corporation nor the names of the
// // contributors may be used to endorse or promote products derived from
// // this software without specific prior written permission.
// //
// // THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// // EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// // PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// // CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// // PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// // LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// // NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// //
// // Questions? Contact
// //                    Jonathan Hu       (jhu@sandia.gov)
// //                    Andrey Prokopenko (aprokop@sandia.gov)
// //                    Ray Tuminaro      (rstumin@sandia.gov)
// //
// // ***********************************************************************
// //
// // @HEADER
/*
 * Xpetra_RegionAMG_decl.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */

#ifndef XPETRA_REGIONAMG_DECL_HPP
#define XPETRA_REGIONAMG_DECL_HPP

// Xpetra
#include <Xpetra_Operator.hpp>
#include "Xpetra_MatrixSplitting.hpp"
#include "Xpetra_Level_def.hpp"

// MueLu
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

namespace Xpetra {

/*!
 * class RegionAMG
 *
 * This class inherits from Xpetra::Operator because the goal is to allow its usage
 * as a preconditioner for a Belos-type linear solver
 *
 * RegionAMG is a class used to split a composite stiffness matrix into region matrices.
 * Each region matrix is associated with a geometric region of the discretized domain.
 * Therefore, the splitting operates in a domain decomposition fashion.
 * Since there are infinitely many ways to split a given composite matrix,
 * we explains here below the rationale behind the approach adopted in this class.
 *  */

template <class Scalar        = MultiVector<>::scalar_type,
          class LocalOrdinal  = typename MultiVector<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class RegionAMG : Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal> multivector_type;
  typedef Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef MatrixSplitting<Scalar, LocalOrdinal, GlobalOrdinal, Node, Xpetra::UseTpetra, false> tpetra_splitting;
  typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> Hierarchy;
  typedef MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> mv_factory_type;
  typedef Level<Scalar, LocalOrdinal, GlobalOrdinal, Node> level;

 public:
  //! Constructors
  //@{

  RegionAMG() {
    std::cout << "This version of constructor is not implemented yet" << std::endl;
  };

  //! This is the actual constructor we are interested in calling
  RegionAMG(Teuchos::RCP<tpetra_splitting> matrixSplitting,
            Teuchos::RCP<Xpetra::RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionHandler,
            RCP<const Teuchos::Comm<int> > comm, Teuchos::ParameterList muelu,
            GlobalOrdinal num_levels, GlobalOrdinal coarsening_factor);

  //@}

  GlobalOrdinal GetNumLevels() { return num_levels_; }

  //! Methods to extract Map information
  //@{
  //! For now, the domain Map coincides with the Domain Map of the composite matrix at the fine level
  virtual RCP<const map_type> getDomainMap() const { return domainMap_; }

  //! For now, the domain Map coincides with the Range Map of the composite matrix at the fine level
  virtual RCP<const map_type> getRangeMap() const { return rangeMap_; }

  //@}

  //! Apply method
  //@{
  //! N.B.: The implementation still has to be finished
  virtual void
  apply(const multivector_type& X, multivector_type& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;
  //@}

  //! hasTransposeApply should not be needed
  virtual bool hasTransposeApply() const { return false; }

 private:
  //! Private variables
  //@{

  // Total number of levels in the hierarchy
  GlobalOrdinal num_levels_ = -1;

  // Total number of regions
  GlobalOrdinal num_regions_ = -1;

  // Coarsening factor to transfer quantities across levels
  int coarsening_factor_ = -1;

  RCP<const Teuchos::Comm<int> > comm_;

  Teuchos::ParameterList muelu_;

  RCP<const map_type> domainMap_;
  RCP<const map_type> rangeMap_;

  // matrixSplitting associated with the composite matrix at the fine level
  RCP<tpetra_splitting> matrixSplitting_;

  // Array of MueLu hierarchies (one for each region)
  Array<RCP<Hierarchy> > regionHierarchies_;

  // Array of levels in the new hierarchy (each levels contains quantities associaed with every region for that level)
  Array<RCP<level> > levels_;

  //@}

  //@Private Methods
  //
  // This methods construct the Hierarchy using the quantities stored in the MueLu::Hierarchy objects
  // There are as many MueLu::Hierarchy objects as the number of geometric regions partitioning the domain
  void SetUpHierarchy();

  // This method is called in SetUpHierarchy and it is in chanrge of extracting quantities from MueLu::Hierarchy objects
  // and store them inside Xpetra::Level objects
  void DefineLevels();

  // Method to create region input multivectors from the composite input multivector
  void computeRegionX(const multivector_type&, Array<RCP<multivector_type> >) const;

  // Method to create composite output multivector from region output multivectors
  void computeCompositeY(Array<RCP<const multivector_type> >, multivector_type&) const;

  // This method detects entries in region output mutlivectors associated with mesh nodes that lie on interregion interfaces
  // and it scales them. The scaling factor is equal to the number of regions that the given degrees of freedom are shared by
  void rescaleInterfaceEntries(Array<RCP<multivector_type> >) const;

  virtual void regionToAllCoarsen(const level&, level&);

  //@}
};

}  // namespace Xpetra

#endif
