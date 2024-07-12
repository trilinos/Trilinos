// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Xpetra_Level_def.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_LEVEL_DEF_HPP
#define XPETRA_LEVEL_DEF_HPP

#include "Xpetra_Level_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Level(int levelID, int num_regions)
  : levelID_(levelID)
  , num_regions_(num_regions) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetP(Array<RCP<matrix_type> >& P) {
  TEUCHOS_TEST_FOR_EXCEPTION(regionP_.size() != 0, Exceptions::RuntimeError, "Current level already has prolongators \n");
  TEUCHOS_TEST_FOR_EXCEPTION(P.size() != num_regions_, Exceptions::RuntimeError, "Number of region prolongators is " << P.size() << "and does not math the number of regions which is " << num_regions_ << " \n");
  regionP_ = P;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetR(Array<RCP<matrix_type> >& R) {
  TEUCHOS_TEST_FOR_EXCEPTION(regionR_.size() != 0, Exceptions::RuntimeError, "Current level already has restrictions \n");
  TEUCHOS_TEST_FOR_EXCEPTION(R.size() != num_regions_, Exceptions::RuntimeError, "Number of region restrictions is " << R.size() << "and does not math the number of regions which is " << num_regions_ << " \n");
  regionR_ = R;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetA(Array<RCP<matrix_type> >& A) {
  TEUCHOS_TEST_FOR_EXCEPTION(regionA_.size() != 0, Exceptions::RuntimeError, "Current level already has operators \n");
  TEUCHOS_TEST_FOR_EXCEPTION(A.size() != num_regions_, Exceptions::RuntimeError, "Number of region operators is " << A.size() << "and does not math the number of regions which is " << num_regions_ << " \n");
  regionA_ = A;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSmoother(Array<RCP<multivector_type> >& S) {
  TEUCHOS_TEST_FOR_EXCEPTION(regionSmoother_.size() != 0, Exceptions::RuntimeError, "Current level already has smoothers \n");
  TEUCHOS_TEST_FOR_EXCEPTION(S.size() != num_regions_, Exceptions::RuntimeError, "Number of region smoothers is " << S.size() << "and does not math the number of regions which is " << num_regions_ << " \n");
  regionSmoother_ = S;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetRegionToAll(Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > regionToAll) {
  TEUCHOS_TEST_FOR_EXCEPTION(regionToAll.size() != num_regions_, Exceptions::RuntimeError, "Passed regionToAll has number of regions equal to " << regionToAll.size() << "which does not match the number of regions in level ID: <<" << levelID_ << " with number of regions equal to " << num_regions_ << "\n");
  level_regionToAll_ = regionToAll;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetNumRegions() const {
  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ == 0, Exceptions::RuntimeError, "level ID: <<" << levelID_ << " does NOT have defined regions \n");
  return num_regions_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetRegionToAll() const {
  TEUCHOS_TEST_FOR_EXCEPTION(level_regionToAll_.size() == 0, Exceptions::RuntimeError, "level ID: " << levelID_ << " does NOT have level_regionToAll_ initialized yet \n");
  return level_regionToAll_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetCompositeIndex(int region_idx, GlobalOrdinal region_node_idx) const {
  TEUCHOS_TEST_FOR_EXCEPTION(level_regionToAll_.size() == 0, Exceptions::RuntimeError, "level ID: " << levelID_ << " does NOT have level_regionToAll_ initialized yet \n");
  TEUCHOS_TEST_FOR_EXCEPTION(level_regionToAll_.size() != num_regions_, Exceptions::RuntimeError, "level ID: " << levelID_ << " has information stored for a number of regions that does NOT match with the declared number of regions \n");
  TEUCHOS_TEST_FOR_EXCEPTION(region_idx >= num_regions_, Exceptions::RuntimeError, "level ID: " << levelID_ << " Invalid region index \n");

  GlobalOrdinal composite_index = -1;

  checkerRegionToAll<GlobalOrdinal> unaryPredicate(region_node_idx);
  typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator global_iterator;
  Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > regionToAll = level_regionToAll_[region_idx];
  global_iterator                                              = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerRegionToAll<GlobalOrdinal> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
  TEUCHOS_TEST_FOR_EXCEPTION(global_iterator == level_regionToAll_[region_idx].end(), Exceptions::RuntimeError, " - Region: " << region_idx << " - "
                                                                                                                              << " node with region index: " << region_idx << " is not in regionToAll[" << region_idx << "]"
                                                                                                                              << "\n");
  composite_index = std::get<1>(*global_iterator);

  return composite_index;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkConsistency() const {
  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ <= 0, Exceptions::RuntimeError, "level ID: " << levelID_ << " does not have any regions \n");
  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ != regionA_.size(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - number of region matrices " << regionA_.size() << " does NOT match the number of regions " << num_regions_ << "\n");
  TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ != level_regionToAll_.size(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - size of level_regionToAll_ " << level_regionToAll_.size() << " does NOT match the number of regions " << num_regions_ << "\n");

  if (levelID_ > 0) {
    TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ != regionP_.size(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - number of prolongators " << regionP_.size() << " does NOT match the number of regions " << num_regions_ << "\n");
    TEUCHOS_TEST_FOR_EXCEPTION(num_regions_ != regionR_.size(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - number of restrictions " << regionR_.size() << " does NOT match the number of regions " << num_regions_ << "\n");
  }

  for (int region_idx = 0; region_idx < num_regions_; ++region_idx) {
    if (levelID_ > 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(regionP_[region_idx]->getGlobalNumCols() != regionA_[region_idx]->getGlobalNumCols(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - Region: " << region_idx << " , numCols(P) = " << regionP_[region_idx]->getGlobalNumCols() << " BUT numCols(A) = " << regionA_[region_idx]->getGlobalNumCols() << "\n");
      TEUCHOS_TEST_FOR_EXCEPTION(regionR_[region_idx]->getGlobalNumRows() != regionA_[region_idx]->getGlobalNumCols(), Exceptions::RuntimeError, "level ID: " << levelID_ << " - Region: " << region_idx << " , numRows(R) = " << regionR_[region_idx]->getGlobalNumRows() << " BUT numCols(A) = " << regionA_[region_idx]->getGlobalNumCols() << "\n");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(regionA_[region_idx]->getRowMap()->getLocalNumElements() > 0 && level_regionToAll_[region_idx].size() != regionA_[region_idx]->getGlobalNumRows(), Exceptions::RuntimeError, "Process ID: " << regionA_[region_idx]->getRowMap()->getComm()->getRank() << " - level ID: " << levelID_ << " - Region: " << region_idx << " , size(regionToAll) = " << level_regionToAll_[region_idx].size() << " BUT numCols(A) = " << regionA_[region_idx]->getGlobalNumCols() << "\n");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Level<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeRegionJacobi() {
  regionSmoother_.resize(num_regions_);

  for (int i = 0; i < num_regions_; ++i) {
    regionSmoother_[i] = VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(regionA_[i]->getRowMap());
    regionA_[i]->getLocalDiagCopy(*(regionSmoother_[i]));
  }
}

}  // namespace Xpetra
#endif
