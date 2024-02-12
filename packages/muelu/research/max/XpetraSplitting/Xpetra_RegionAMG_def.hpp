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
 * Xpetra_RegionAMG_def.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_REGIONAMG_DEF_HPP
#define XPETRA_REGIONAMG_DEF_HPP

#include "Xpetra_RegionAMG_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RegionAMG(
    Teuchos::RCP<tpetra_splitting> matrixSplitting,
    Teuchos::RCP<Xpetra::RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionHandler,
    RCP<const Teuchos::Comm<int> > comm, Teuchos::ParameterList muelu,
    GlobalOrdinal num_levels, GlobalOrdinal coarsening_factor)
  : num_levels_(num_levels)
  , coarsening_factor_(coarsening_factor)
  , comm_(comm)
  , muelu_(
        muelu) {
  TEUCHOS_TEST_FOR_EXCEPT(matrixSplitting.is_null());
  matrixSplitting_ = matrixSplitting;

  // The maps defined here below are used to interface with outer applications (e.g. Belos-type solvers)
  // and guarantee that DomainMap and RangeMap partition input multivectors and output multivectors in a proper way.
  // the maps stored in doaminMap and rangeMap coincide with the maps of the composite matrix
  domainMap_   = matrixSplitting_->getDomainMap();
  rangeMap_    = matrixSplitting_->getRangeMap();
  num_regions_ = matrixSplitting_->getNumRegions();

  SetUpHierarchy();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetUpHierarchy() {
  TEUCHOS_TEST_FOR_EXCEPTION(num_levels_ <= 0, Exceptions::RuntimeError, "Number of levels must be a positive integer number \n");
  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ <= 0, Exceptions::RuntimeError, "No existing regions \n");

  // Creation of the MueLu list for the region multigrid preconditioner
  // It would be ideal to have this parameter list read by an .xml file
  // at the final stage of the implementation.
  // As for now, we hard code it insice of SetUpHierarchy because we want standard setups for almost everything
  RCP<ParameterList> list = rcp(new Teuchos::ParameterList());
  list->setName("MueLu");

  list->set("verbosity", "low");
  list->set("number of equations", 1);
  list->set("max levels", num_levels_);
  list->set("multigrid algorithm", "unsmoothed");
  list->set("smoother: pre or post", "none");
  list->set("coarse: type", "RELAXATION");
  list->set("coarse: max size", 16);

  GlobalOrdinal num_regions = matrixSplitting_->getNumRegions();
  regionHierarchies_.resize(num_regions);

  // Instantiation of MueLu Hierarchies for each region of the domain
  for (GlobalOrdinal region_idx = 0; region_idx < num_regions; ++region_idx) {
    // The parameter n is the total number of mesh nodes in a region
    GlobalOrdinal n = matrixSplitting_->getRegionHandler()->GetNumRegionNodes(region_idx);

    // N.B. The following instructions is based on the assumption that each region is a square. If not, then another instruciton
    // must be used to specify the size of the region mesh along teh x-direction
    GlobalOrdinal nx = std::sqrt(n);

    RCP<const map_type> map                                  = matrixSplitting_->getRegionMatrix(region_idx)->getRowMap();
    size_t NumMyElements                                     = map->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

    RCP<multivector_type> coords = mv_factory_type::Build(map, 2);

    // Coordinates are created to pass them as additional input the MueLu::Hierarchy constructors
    // The only thing we care aboput is to create a unique mapping between coordinates and degrees of freedom associated with mesh nodes
    // Therefore we base the creation of coordinates upon a fake region mesh, which does not coincide with its real gemetric shape
    // The fake region we use for coordinates is a [0,(number_nodes_x-1)]X[0,(number_nodes_y-1)] rectangle
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > Coord(2);
    Coord[0] = coords->getDataNonConst(0);
    Coord[1] = coords->getDataNonConst(1);

    for (size_t i = 0; i < NumMyElements; ++i) {
      GlobalOrdinal ix = MyGlobalElements[i] % nx;
      GlobalOrdinal iy = (MyGlobalElements[i] - ix) / nx;

      Coord[0][i] = Teuchos::as<Scalar>(ix);
      Coord[1][i] = Teuchos::as<Scalar>(iy);
    }

    // We create as many MueLu::Hierarchy objects as the nubmer of region that partition the composite domain
    regionHierarchies_[region_idx] = MueLu::CreateXpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(matrixSplitting_->getRegionMatrix(region_idx), *list, coords);

    // // The following commented line here below should replace the line right above if the MueLu parameter list is passed through an .xml file
    // regionHierarchies_[region_idx] = MueLu::CreateXpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>( matrixSplitting_->getRegionMatrix( region_idx), muelu_, coords );

    // Different regions may have meshes with different size, or the number of levels the user wants to create may be too big
    // with respect to the number of levels MueLu allows (i.e. the minimum coarse size may be reached for a smaller number of levels)
    // In this case, the value of num_levels_ is readjusted to the minimum number of levels instantied across all the regions of the domain
    num_levels_ = std::min(num_levels_, regionHierarchies_[region_idx]->GetNumLevels());
  }

  DefineLevels();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DefineLevels() {
  levels_.clear();
  Array<RCP<matrix_type> > P;
  P.clear();
  Array<RCP<matrix_type> > R;
  R.clear();
  Array<RCP<matrix_type> > A;
  A.clear();

  TEUCHOS_TEST_FOR_EXCEPTION(levels_.size() != 0, Exceptions::RuntimeError, "Levels structure is already initialized \n");
  for (int i = 0; i < num_levels_; ++i) {
    P.clear();
    R.clear();
    A.clear();

    RCP<level> new_level = rcp(new level(i, num_regions_));
    for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
      TEUCHOS_TEST_FOR_EXCEPTION(!regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("A"), Exceptions::RuntimeError, "No existing operator at level " << i << " of region " << region_idx << "\n");
      A.push_back(regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("A"));

      // MueLu::Hierarchy objects store prolongator and restriction operators only in levels that are the result of a coarsening (they do not exist at the fine level)
      if (i > 0) {
        TEUCHOS_TEST_FOR_EXCEPTION(!regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("P"), Exceptions::RuntimeError, "No existing prolongator at level " << i << " of region " << region_idx << "\n");
        TEUCHOS_TEST_FOR_EXCEPTION(!regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("R"), Exceptions::RuntimeError, "No existing restriction at level " << i << " of region " << region_idx << "\n");
        P.push_back(regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("P"));
        R.push_back(regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("R"));
      }
    }
    // Following the same policy adopted in MueLu::Hierarchy, prolongator and restriction operators are stored at levels produced with a coarsening
    if (i > 0) {
      new_level->SetP(P);
      new_level->SetR(R);
    }

    // The Galerkin operator is define dat every level (regardless of coarsening)
    new_level->SetA(A);

    // The structure regionToAll is passed to the fine level as it is
    if (0 == i)
      new_level->SetRegionToAll(matrixSplitting_->getRegionHandler()->GetRegionToAll());
    else {
      // Coarsening level own a regionToAll structure whcih is the result of an injected coarsened regionToAll from the upper level of the hierarchy
      Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > coarse_regionToAll;
      regionToAllCoarsen(*(levels_[i - 1]), *new_level);
    }

    new_level->checkConsistency();

    // // The following routine should be dedicated to the construction of the region smoother (where it could be possible to define
    // the region smoother as a local visualization of the composite smoother)
    // FIXME ComputeRegionJacobi AS FOR NOW LOCAL JACOBI IS THE ACTUAL REGIONAL SMOOTHER
    // (THIS HAS TO BE CHANGED SO THAT A PORTION OF THE COMPOSITE SMOOTHER IS STORED INSTEAD)
    new_level->ComputeRegionJacobi();

    levels_.push_back(new_level);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::regionToAllCoarsen(const level& fine, level& coarse) {
  TEUCHOS_TEST_FOR_EXCEPTION(fine.GetNumRegions() != coarse.GetNumRegions(), Exceptions::RuntimeError, "Level " << fine.GetLevelID() << "has " << fine.GetNumRegions() << " regions instantiated whereas level " << coarse.GetLevelID() << "has " << coarse.GetNumRegions() << " regions instantiated \n");

  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > fine_regionToAll = fine.GetRegionToAll();
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > coarse_regionToAll(num_regions_);

  for (GlobalOrdinal region_idx = 0; region_idx < num_regions_; ++region_idx) {
    // Coordinates needed to apply the coarsening
    GlobalOrdinal n  = fine.GetNumRegionNodes(region_idx);
    GlobalOrdinal nx = std::sqrt(n);
    GlobalOrdinal ny = nx;

    TEUCHOS_TEST_FOR_EXCEPTION((nx - 1) % coarsening_factor_ != 0 || (ny - 1) % coarsening_factor_ != 0, Exceptions::RuntimeError, "Region: " << region_idx << " cannot be coarsened by factor equal to " << coarsening_factor_ << "n: " << n << " - nx: " << nx << " - ny: " << ny << " \n");

    // here below we exploit the fact that the regionToAll has been sorted in ascending order for the region node index
    GlobalOrdinal count = 1;

    while (count <= fine_regionToAll[region_idx].size()) {
      coarse_regionToAll[region_idx].push_back(fine_regionToAll[region_idx][count]);
      if (count % ny != 0)
        count = count + coarsening_factor_;
      else
        count = count + (coarsening_factor_ - 1) * ny + 1;
    }

    coarse.SetRegionToAll(coarse_regionToAll);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::computeRegionX(const multivector_type& X, Array<RCP<multivector_type> > regionX) const {
  // Array to store extended region Maps (needed to copy composite entries of the input vector into region partitioning of it)
  Teuchos::Array<Teuchos::Array<GlobalOrdinal> > overlapping_composite_array;

  // Create Overlapping composite maps
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    Teuchos::Array<GlobalOrdinal> overlapped_composite;
    LocalOrdinal num_elements = regionX[region_idx]->getMap()->getLocalNumElements();
    for (LocalOrdinal local_region_index = 0; local_region_index < num_elements; ++local_region_index) {
      GlobalOrdinal composite_region_index    = regionX[region_idx]->getMap()->getGlobalElement(local_region_index);
      GlobalOrdinal composite_composite_index = levels_[0]->GetCompositeIndex(region_idx, composite_region_index + 1);
      overlapped_composite.push_back(composite_composite_index - 1);
    }
    overlapping_composite_array.push_back(overlapped_composite);
  }

  // We split at first the input and output multivectors into region ones
  Array<RCP<multivector_type> > composite_overlapping_X;
  composite_overlapping_X.resize(num_regions_);
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    // The new Map associated with a process must contain entries currently owned plus new ones previously owned only by neighbouring processes
    Teuchos::Array<GlobalOrdinal> aux;
    if (overlapping_composite_array[region_idx].size() > 0) {
      Teuchos::Array<GlobalOrdinal> aux1 = overlapping_composite_array[region_idx];
      Teuchos::Array<GlobalOrdinal> aux2 = X.getMap()->getLocalElementList();
      std::set_union(aux1.begin(), aux1.end(), aux2.begin(), aux2.end(), std::back_inserter(aux));
    } else {
      aux = X.getMap()->getLocalElementList();
    }
    RCP<map_type> overlapping_composite_map                 = MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Xpetra::UseTpetra, X.getGlobalLength(), aux, 0, comm_);
    composite_overlapping_X[region_idx]                     = MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(overlapping_composite_map, X.getNumVectors());
    RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Import1 = ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(X.getMap(), overlapping_composite_map);
    TEUCHOS_TEST_FOR_EXCEPTION(X.getMap()->getMinAllGlobalIndex() != composite_overlapping_X[region_idx]->getMap()->getMinAllGlobalIndex(), Exceptions::RuntimeError, "Minimal index in old an new maps do not coincide \n");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getMap()->getMaxAllGlobalIndex() != composite_overlapping_X[region_idx]->getMap()->getMaxAllGlobalIndex(), Exceptions::RuntimeError, "Maximal index in old an new maps do not coincide \n");
    composite_overlapping_X[region_idx]->doImport(X, *Import1, Xpetra::INSERT);
  }

  // Copy values from composite input multivector X into region multivectors regionX
  for (int i = 0; i < X.getNumVectors(); ++i) {
    for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
      ArrayRCP<const Scalar> composite_column = composite_overlapping_X[region_idx]->getData(i);
      LocalOrdinal num_elements               = regionX[region_idx]->getMap()->getLocalNumElements();
      ArrayRCP<Scalar> region_column          = regionX[region_idx]->getDataNonConst(i);
      for (LocalOrdinal local_region_index = 0; local_region_index < num_elements; ++local_region_index) {
        GlobalOrdinal composite_region_index    = regionX[region_idx]->getMap()->getGlobalElement(local_region_index);
        GlobalOrdinal composite_composite_index = levels_[0]->GetCompositeIndex(region_idx, composite_region_index + 1);
        LocalOrdinal local_composite_index      = composite_overlapping_X[region_idx]->getMap()->getLocalElement(composite_composite_index - 1);
        region_column[local_region_index]       = composite_column[local_composite_index];
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::computeCompositeY(Array<RCP<const multivector_type> > regionY, multivector_type& Y) const {
  // Array to store extended region Maps (needed to copy composite entries of the input vector into region partitioning of it)
  Teuchos::Array<GlobalOrdinal> overlapping_composite;

  // Create Overlapping composite maps
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    LocalOrdinal num_elements = regionY[region_idx]->getMap()->getLocalNumElements();
    for (LocalOrdinal local_region_index = 0; local_region_index < num_elements; ++local_region_index) {
      GlobalOrdinal composite_region_index    = regionY[region_idx]->getMap()->getGlobalElement(local_region_index);
      GlobalOrdinal composite_composite_index = levels_[0]->GetCompositeIndex(region_idx, composite_region_index + 1);
      overlapping_composite.push_back(composite_composite_index - 1);
    }
  }

  typename Teuchos::Array<GlobalOrdinal>::iterator last;
  std::sort(overlapping_composite.begin(), overlapping_composite.end());
  last = std::unique(overlapping_composite.begin(), overlapping_composite.end());
  overlapping_composite.erase(last, overlapping_composite.end());

  // We split at first the input and output multivectors into region ones
  RCP<multivector_type> composite_overlapping_Y;
  RCP<map_type> overlapping_composite_map = MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Xpetra::UseTpetra, Y.getGlobalLength(), overlapping_composite, 0, comm_);
  composite_overlapping_Y                 = MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(overlapping_composite_map, Y.getNumVectors());

  // Safety checks to guarantee the legitimacy of non-unique source composite map and unique target composite map
  // The unique target composite map must be a subset of the non-unique source composite map
  // Before running the std::includes algorithm it is very important to sort the ElementList in ascending order
  Teuchos::Array<GlobalOrdinal> original_map    = Y.getMap()->getLocalElementList();
  Teuchos::Array<GlobalOrdinal> overlapping_map = composite_overlapping_Y->getMap()->getLocalElementList();
  std::sort(original_map.begin(), original_map.end());
  std::sort(overlapping_map.begin(), overlapping_map.end());
  TEUCHOS_TEST_FOR_EXCEPTION(!(std::includes(overlapping_map.begin(), overlapping_map.end(), original_map.begin(), original_map.end())), Exceptions::RuntimeError, "Overlapping (non-unique) composite map does NOT include original (unique) composite map \n");

  // Copy values from output region multivectors regionY into output composite multivector Y
  for (int i = 0; i < Y.getNumVectors(); ++i) {
    for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
      ArrayRCP<Scalar> composite_column    = composite_overlapping_Y->getDataNonConst(i);
      LocalOrdinal num_elements            = regionY[region_idx]->getMap()->getLocalNumElements();
      ArrayRCP<const Scalar> region_column = regionY[region_idx]->getData(i);
      for (LocalOrdinal local_region_index = 0; local_region_index < num_elements; ++local_region_index) {
        GlobalOrdinal composite_region_index    = regionY[region_idx]->getMap()->getGlobalElement(local_region_index);
        GlobalOrdinal composite_composite_index = levels_[0]->GetCompositeIndex(region_idx, composite_region_index + 1);
        LocalOrdinal local_composite_index      = composite_overlapping_Y->getMap()->getLocalElement(composite_composite_index - 1);
        composite_column[local_composite_index] += region_column[local_region_index];
      }
    }
  }

  // Import object needed to transfer output composite multivector from non-unique map to unique map
  RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Import1 = ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(composite_overlapping_Y->getMap(), Y.getMap());

  // IMPORTANT: ALWAYS USE IMPORT-EXPORT OBJECTS IN FORWARD MODE (i.e. doImport methods must use Import objects and doExport methods must use Export objects)
  // Reverse mode (i.e. crossed use of doImport/doExport methods with Export/Import objects makes the code possibly crash due to undefined behavior)
  Y.doImport(*composite_overlapping_Y, *Import1, Xpetra::ADD);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::rescaleInterfaceEntries(Array<RCP<multivector_type> > regionY) const {
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > regionToAll    = levels_[0]->GetRegionToAll();
  Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > > interfaceNodes = matrixSplitting_->getRegionHandler()->GetInterfaceNodes();

  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ != regionToAll.size(), Exceptions::RuntimeError, "Regions stored in Level 0 do not match with total number of regions in RegionAMG class \n");

  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    LocalOrdinal num_elements = regionY[region_idx]->getMap()->getLocalNumElements();
    for (LocalOrdinal local_region_index = 0; local_region_index < num_elements; ++local_region_index) {
      GlobalOrdinal composite_region_index    = regionY[region_idx]->getMap()->getGlobalElement(local_region_index);
      GlobalOrdinal composite_composite_index = levels_[0]->GetCompositeIndex(region_idx, composite_region_index + 1);
      checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(composite_composite_index);
      typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
      nodes_to_region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicateNode);
      // The rescaling must be applied only to entries associated with mesh on the interface
      if (nodes_to_region_iterator != interfaceNodes.end()) {
        Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);
        for (int i = 0; i < regionY[region_idx]->getNumVectors(); ++i) {
          ArrayRCP<Scalar> region_column = regionY[region_idx]->getDataNonConst(i);
          if (nodal_regions.size() > 1)
            region_column[local_region_index] = region_column[local_region_index] / (nodal_regions.size());
        }
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const multivector_type& X, multivector_type& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  // N.B.: currently Scalar quantities alpha and beta are passed as input parameters to have the apply method signature match with the apply signature of an Xpetra::Operator
  // however we are not currently using them (for us alpha=0 and beta=1)

  // At first we check that input and output vector have matching maps with the composite matrix
  // The Map of X must coincide with the Domain map of compositeA
  // The Map of Y must coincide with the Range map of compositeA
  TEUCHOS_TEST_FOR_EXCEPTION(!(X.getMap()->isSameAs(*(matrixSplitting_->getMatrix()->getDomainMap()))), Exceptions::RuntimeError, "Map of composite input multivector X does not coincide with Domain Map of composite matrix \n");
  TEUCHOS_TEST_FOR_EXCEPTION(!(Y.getMap()->isSameAs(*(matrixSplitting_->getMatrix()->getRangeMap()))), Exceptions::RuntimeError, "Map of composite input multivector X does not coincide with Range Map of composite matrix \n");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), Exceptions::RuntimeError, "Number of vectors in input numltivector X does NOT match number of vectors in output multivector Y \n");

  // We split at first the input and output multivectors into region ones
  Array<RCP<multivector_type> > regionX;
  Array<RCP<multivector_type> > regionY;
  regionX.resize(num_regions_);
  regionY.resize(num_regions_);

  // Associate Maps to region input (regionX) and output vectors (regionY)
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    regionX[region_idx] = MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(levels_[0]->GetRegionMatrix(region_idx)->getDomainMap(), X.getNumVectors());
    regionY[region_idx] = MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(levels_[0]->GetRegionMatrix(region_idx)->getRangeMap(), Y.getNumVectors());
  }

  // Split the composite input multivector	X into region multivector regionX
  computeRegionX(X, regionX);

  // This is the portion where the V-cycle is executed (for now we only have region V-cycles)
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    regionHierarchies_[region_idx]->Iterate(*regionX[region_idx], *regionY[region_idx]);
  }

  // We rescale entries of region multivector regionY that are associated with mesh nodes on an interface
  rescaleInterfaceEntries(regionY);

  // We create const view of the region multivector regionY because we want to guarantee that the composite output multivector Y
  // is constructed without modifying any information given from each region
  Array<RCP<const multivector_type> > regionYconst;
  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    regionYconst.push_back(regionY[region_idx].getConst());
  }

  // Assemble composite output multivector Y from region multivector regionY
  computeCompositeY(regionYconst, Y);
}

}  // namespace Xpetra
#endif
