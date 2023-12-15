// @HEADER
//
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_Q2Q1UPFACTORY_DECL_HPP
#define MUELU_Q2Q1UPFACTORY_DECL_HPP

#include <Teuchos_RCP.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"

#include <algorithm>
#include <vector>
#include <sys/stat.h>

namespace MueLu {

template <class LocalOrdinal = int>
class MyCptList_ {
  typedef LocalOrdinal LO;

 public:
  MyCptList_(int n, int nnzPerRow = 100) {
    TEUCHOS_TEST_FOR_EXCEPTION(nnzPerRow <= 0, Exceptions::RuntimeError, "Why is nnzPerRow " << nnzPerRow << "?");

    nnzPerRow_ = nnzPerRow;
    storage_.resize(n * nnzPerRow);
    numCpts_.resize(n, 0);

    list_.resize(n, NULL);
    for (int i = 0; i < n; i++)
      list_[i] = &storage_[i * nnzPerRow];
  }

  size_t getLocalNumRows() const { return list_.size(); }
  int getNnzPerRow() const { return nnzPerRow_; }
  std::vector<short>& getNumCpts() { return numCpts_; }
  Teuchos::Array<LO>& getCList() { return cptlist_; }
  const Teuchos::Array<LO>& getCList() const { return cptlist_; }
  const std::vector<short>& getNumCpts() const { return numCpts_; }
  LO* operator()(int i) { return list_[i]; }
  const LO* operator()(int i) const { return list_[i]; }

 private:
  std::vector<LO*> list_;       // list[k] gives the CPOINTs that interpolate to the k-th fine point
                                // These CPOINTs are given as fine grid *local* indices
  std::vector<short> numCpts_;  // Number of CPOINTs for each point
  Teuchos::Array<LO> cptlist_;  // List of CPOINTs for each point
  int nnzPerRow_;               // Max number of CPOINTs per row in order for a row to use storage_
  std::vector<LO> storage_;     // Large data array used to store most CPOINT information
};

template <class Scalar        = MueLu::DefaultScalar,
          class LocalOrdinal  = MueLu::DefaultLocalOrdinal,
          class GlobalOrdinal = MueLu::DefaultGlobalOrdinal,
          class Node          = MueLu::DefaultNode>
class Q2Q1uPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"
  typedef MyCptList_<LocalOrdinal> MyCptList;

 private:
  enum Status {
    UNASSIGNED = '0',
    CANDIDATE  = '1',
    FPOINT     = '2',
    TWOTIMER   = '3',
    CPOINT     = '4',
    CPOINT_U   = '5'
  };

  std::string getStatusString(char status) const {
    switch (status) {
      case UNASSIGNED: return "UNASSIGNED";
      case CANDIDATE: return "CANDIDATE";
      case FPOINT: return "FPOINT";
      case TWOTIMER: return "TWOTIMER";
      case CPOINT: return "CPOINT";
      case CPOINT_U: return "CPOINT_U";
      default: return "UNKNOWN";
    }
  }

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  Q2Q1uPFactory() {}

  //! Destructor.
  virtual ~Q2Q1uPFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

 private:
  void FindDist4Cpts(const Matrix& A, const MultiVector& coords, const Array<LO>& userCpts, std::vector<char>& status, MyCptList& myCpts, int levelID) const;
  void PhaseTwoPattern(const Matrix& A, const MultiVector& coords, const std::vector<char>& status, MyCptList& myCpts) const;
  void FindMidPoints(const Matrix& A, const MultiVector& coords, Array<LO>& Cptlist, const MyCptList& myCpts) const;
  void CompDistances(const Matrix& A, LO start, int numDist, std::vector<LO>& dist1, std::vector<LO>& dist2,
                     std::vector<LO>& dist3, std::vector<LO>& dist4) const;
  void CreateCrsPointers(const Matrix& A, ArrayRCP<const size_t>& ia, ArrayRCP<const LO>& ja) const;
  void CptDepends2Pattern(const Matrix& A, const MyCptList& myCpts, RCP<Matrix>& P, LO offset) const;

  void DumpStatus(const std::string& filename, const std::vector<char>& status, int NDim, bool isAmalgamated = true) const;
  void DumpCoords(const std::string& filename, const MultiVector& coords) const;
};

// Sort a double array and move along list2 to match sorted array
template <class SC>
void Muelu_az_dsort2(std::vector<SC>& dlist, std::vector<int>& list2) {
  int l, r, j, i, flag;
  int RR2;
  SC dRR, dK;

  int N = dlist.size();
  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  dRR = dlist[l - 1];
  dK  = dlist[l - 1];

  if (list2.size()) {
    RR2 = list2[l - 1];
    while (r != 0) {
      j    = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (dlist[j] > dlist[j - 1]) j = j + 1;

          if (dlist[j - 1] > dK) {
            dlist[i - 1] = dlist[j - 1];
            list2[i - 1] = list2[j - 1];
          } else {
            flag = 0;
          }
        }
      }
      dlist[i - 1] = dRR;
      list2[i - 1] = RR2;

      if (l == 1) {
        dRR      = dlist[r];
        RR2      = list2[r];
        dK       = dlist[r];
        dlist[r] = dlist[0];
        list2[r] = list2[0];
        r        = r - 1;
      } else {
        l   = l - 1;
        dRR = dlist[l - 1];
        RR2 = list2[l - 1];
        dK  = dlist[l - 1];
      }
    }
    dlist[0] = dRR;
    list2[0] = RR2;
  } else {
    while (r != 0) {
      j    = l;
      flag = 1;
      while (flag == 1) {
        i = j;
        j = j + j;
        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (dlist[j] > dlist[j - 1]) j = j + 1;
          if (dlist[j - 1] > dK) {
            dlist[i - 1] = dlist[j - 1];
          } else {
            flag = 0;
          }
        }
      }
      dlist[i - 1] = dRR;
      if (l == 1) {
        dRR      = dlist[r];
        dK       = dlist[r];
        dlist[r] = dlist[0];
        r        = r - 1;
      } else {
        l   = l - 1;
        dRR = dlist[l - 1];
        dK  = dlist[l - 1];
      }
    }
    dlist[0] = dRR;
  }
}

/* ******************************************************************* */
/* sort an array and move along list2 and/or list to match sorted array*/
/* ------------------------------------------------------------------- */
template <class SC>
void Muelu_az_sort(int list[], int N, int list2[], SC list3[]) {
  int l, r, RR, K, j, i, flag;
  int RR2;
  SC RR3;

  if (N <= 1) return;

  l  = N / 2 + 1;
  r  = N - 1;
  l  = l - 1;
  RR = list[l - 1];
  K  = list[l - 1];

  if ((list2 != NULL) && (list3 != NULL)) {
    RR2 = list2[l - 1];
    RR3 = list3[l - 1];
    while (r != 0) {
      j    = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[i - 1]  = list[j - 1];
            list2[i - 1] = list2[j - 1];
            list3[i - 1] = list3[j - 1];
          } else {
            flag = 0;
          }
        }
      }

      list[i - 1]  = RR;
      list2[i - 1] = RR2;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list[r];
        RR2 = list2[r];
        RR3 = list3[r];

        K        = list[r];
        list[r]  = list[0];
        list2[r] = list2[0];
        list3[r] = list3[0];
        r        = r - 1;
      } else {
        l   = l - 1;
        RR  = list[l - 1];
        RR2 = list2[l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[0]  = RR;
    list2[0] = RR2;
    list3[0] = RR3;
  } else if (list2 != NULL) {
    RR2 = list2[l - 1];
    while (r != 0) {
      j    = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[i - 1]  = list[j - 1];
            list2[i - 1] = list2[j - 1];
          } else {
            flag = 0;
          }
        }
      }

      list[i - 1]  = RR;
      list2[i - 1] = RR2;

      if (l == 1) {
        RR  = list[r];
        RR2 = list2[r];

        K        = list[r];
        list[r]  = list[0];
        list2[r] = list2[0];
        r        = r - 1;
      } else {
        l   = l - 1;
        RR  = list[l - 1];
        RR2 = list2[l - 1];
        K   = list[l - 1];
      }
    }

    list[0]  = RR;
    list2[0] = RR2;
  } else if (list3 != NULL) {
    RR3 = list3[l - 1];
    while (r != 0) {
      j    = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[i - 1]  = list[j - 1];
            list3[i - 1] = list3[j - 1];
          } else {
            flag = 0;
          }
        }
      }

      list[i - 1]  = RR;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list[r];
        RR3 = list3[r];

        K        = list[r];
        list[r]  = list[0];
        list3[r] = list3[0];
        r        = r - 1;
      } else {
        l   = l - 1;
        RR  = list[l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[0]  = RR;
    list3[0] = RR3;

  } else {
    while (r != 0) {
      j    = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[i - 1] = list[j - 1];
          } else {
            flag = 0;
          }
        }
      }

      list[i - 1] = RR;

      if (l == 1) {
        RR = list[r];

        K       = list[r];
        list[r] = list[0];
        r       = r - 1;
      } else {
        l  = l - 1;
        RR = list[l - 1];
        K  = list[l - 1];
      }
    }

    list[0] = RR;
  }
}

// Merge two already sorted lists into one combined sorted list.
// NOTE: lists are given as integer arrays. These integer arrays give
// locations in CoordDist[] defining the list values. That the ith value
// associated with the Candidates list is actually CoordDist[Candidates[i]].
template <class T>
void MergeSort(std::vector<int>& oldCandidates, size_t numOldCandidates, const std::vector<int>& newCandidates, const std::vector<T>& coordDist, ArrayRCP<const size_t> ia) {
  size_t numNewCandidates = newCandidates.size();
  size_t numCandidates    = numOldCandidates + numNewCandidates;

  oldCandidates.resize(numCandidates);

  int i = numOldCandidates - 1;
  int j = numNewCandidates - 1;
  int k = numCandidates - 1;
  while ((i >= 0) || (j >= 0)) {
    if (i < 0)
      oldCandidates[k--] = newCandidates[j--];
    else if (j < 0)
      oldCandidates[k--] = oldCandidates[i--];
    else {
      int ii = oldCandidates[i];
      int jj = newCandidates[j];

      // Must match code above. There is something arbitrary
      // and crappy about the current weighting.

#ifdef optimal
      if (-coordDist[ii] - .01 * (ia[ii + 1] - ia[ii]) + 1.e-10 * (ii + 1) <
          -coordDist[jj] - .01 * (ia[jj + 1] - ia[jj]) + 1.e-10 * (jj + 1))
#else
      if (coordDist[ii] - .0 * (ia[ii + 1] - ia[ii]) + 1.e-3 * (ii + 1) <
          coordDist[jj] - .0 * (ia[jj + 1] - ia[jj]) + 1.e-3 * (jj + 1))
      // if (ii < jj)
#endif
        oldCandidates[k--] = newCandidates[j--];
      else
        oldCandidates[k--] = oldCandidates[i--];
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A");

  RCP<const FactoryBase> rcpThis = rcpFromRef(*this);
  validParamList->set<RCP<const FactoryBase> >("CoordinatesVelocity", rcpThis, "Generating factory of the coordinates");
  validParamList->set<RCP<const FactoryBase> >("AForPat", rcpThis, "Generating factory for Apattern");
  validParamList->set<RCP<const FactoryBase> >("CoordinatesPressure", rcpThis, "Generating factory of the coordinates");
  validParamList->set<RCP<const FactoryBase> >("p2vMap", rcpThis, "Mapping of pressure coords to u-velocity coords");

  validParamList->set<std::string>("mode", "pressure", "Mode");
  validParamList->set<bool>("phase2", false, "Use extra phase to improve pattern");
  validParamList->set<bool>("dump status", false, "Output status");

  validParamList->set<double>("tau_2", sqrt(0.0015), "tau_2 parameter from the paper (used for mid points)");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");

  const ParameterList& pL = GetParameterList();
  bool pressureMode       = (pL.get<std::string>("mode") == "pressure");

  // NOTE: we cannot simply do Input(fineLevel, "CoordinatePressure", as in
  // valid parameter list we specified *this as the generating factory
  if (fineLevel.GetLevelID()) {
    if (pressureMode) {
      Input(fineLevel, "CoordinatesPressure");
    } else {
      Input(fineLevel, "CoordinatesVelocity");
      Input(fineLevel, "AForPat");
      Input(fineLevel, "p2vMap");
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  typedef Teuchos::ScalarTraits<SC> STS;

  const ParameterList& pL = GetParameterList();
  bool pressureMode       = (pL.get<std::string>("mode") == "pressure");
  GetOStream(Runtime0) << (pressureMode ? "Pressure" : "Velocity") << " mode" << std::endl;

  bool fineLevelID = fineLevel.GetLevelID();

  RCP<Matrix> A           = Get<RCP<Matrix> >(fineLevel, "A");
  Xpetra::global_size_t N = A->getRowMap()->getGlobalNumElements();

  RCP<MyCptList> myCpts = rcp(new MyCptList(N));
  std::vector<char> status(N, UNASSIGNED);

  RCP<MultiVector> coords;
  RCP<Matrix> AForPat;
  int NDim = -1;
  if (pressureMode) {
    if (fineLevelID == 0)
      coords = fineLevel.Get<RCP<MultiVector> >("CoordinatesPressure", NoFactory::get());
    else
      coords = Get<RCP<MultiVector> >(fineLevel, "CoordinatesPressure");
    NDim = coords->getNumVectors();

    Array<LO> userCpts;  // pressure does not reuse any CPOINTs
    FindDist4Cpts(*A, *coords, userCpts, status, *myCpts, fineLevelID);

    if (pL.get<bool>("phase2")) {
      // Beef up any limited pattern
      PhaseTwoPattern(*A, *coords, status, *myCpts);
    }

  } else {
    // Do all the coarsening/pattern stuff on amalgamated velocities.
    // We need to guarantee that all velocity dofs are treated identically
    // This means that we must amalgmate AForPat and the velocity coordinates
    if (fineLevelID == 0)
      coords = fineLevel.Get<RCP<MultiVector> >("CoordinatesVelocity", NoFactory::get());
    else
      coords = Get<RCP<MultiVector> >(fineLevel, "CoordinatesVelocity");
    if (fineLevelID == 0)
      AForPat = fineLevel.Get<RCP<Matrix> >("AForPat", NoFactory::get());
    else
      AForPat = Get<RCP<Matrix> >(fineLevel, "AForPat");
    NDim = coords->getNumVectors();

    TEUCHOS_TEST_FOR_EXCEPTION(!coarseLevel.IsAvailable("PresCptsAndMids"), Exceptions::RuntimeError,
                               "Pressure points are not available");

    Array<LO> userCpts = coarseLevel.Get<Array<LO> >("PresCptsAndMids");
    GetOStream(Runtime1) << "Found stored pressure C-points: " << userCpts.size() << " " << userCpts << std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION(N % NDim, Exceptions::RuntimeError, "Number of velocity DOFs is odd");
    Xpetra::global_size_t NN = N / NDim;

    std::vector<GO> gNodeIds(NN);
    for (size_t k = 0; k < NN; k++)
      gNodeIds[k] = k;

    RCP<Map> nodeMap = MapFactory::Build(AForPat->getRowMap()->lib(), NN, gNodeIds, 0, AForPat->getRowMap()->getComm());

    // FIXME: remove magic number 30
    RCP<Matrix> amalgA        = MatrixFactory::Build(nodeMap, nodeMap, 30);
    RCP<CrsMatrix> amalgA_crs = rcp_dynamic_cast<CrsMatrixWrap>(amalgA)->getCrsMatrix();

    // FIXME: this should be written similar to CoalesceDropFactory Merge
    for (LO row = 0; row < as<LO>(AForPat->getRowMap()->getLocalNumElements()); row += NDim) {
      GO grid      = AForPat->getRowMap()->getGlobalElement(row);
      GO currentId = grid / NDim;

      Teuchos::ArrayView<const LO> inds;
      Teuchos::ArrayView<const SC> vals;
      AForPat->getLocalRowView(row, inds, vals);

      size_t nnz = inds.size();

      // Count the number of nonzero block columns in this row
      // NOTE: this assumes that blocks are dense, i.e. that if one column is
      // nonzero, then all columns in the same block are nonzeros
      LO realnnz = 0;
      for (LO col = 0; col < Teuchos::as<LO>(nnz); col++)
        if (inds[col] % NDim == 0)
          realnnz++;

      if (realnnz == 0)
        continue;

      Teuchos::Array<GO> cnodeIds(realnnz, 0);
      Teuchos::Array<SC> ones(realnnz, STS::one());  // Pattern has all 1's

      realnnz = 0;
      for (LO col = 0; col < Teuchos::as<LO>(nnz); col++) {
        if (inds[col] % NDim == 0) {
          GO gcid             = AForPat->getColMap()->getGlobalElement(inds[col]);
          cnodeIds[realnnz++] = gcid / NDim;
        }
      }
      amalgA_crs->insertGlobalValues(currentId, cnodeIds, ones);
    }
    amalgA_crs->fillComplete(nodeMap, nodeMap);

    // Amalgmate the velocity coordinates
    // NOTE: This assumes that the original coords vector contains duplicated (x NDim) degrees of freedom
    RCP<MultiVector> amalgCoords = Xpetra::MultiVectorFactory<SC, LO, GO, Node>::Build(nodeMap, NDim);

    for (int j = 0; j < NDim; j++) {
      ArrayRCP<SC> coordView      = coords->getDataNonConst(j);
      ArrayRCP<SC> amalgCoordView = amalgCoords->getDataNonConst(j);
      for (size_t k = 0; k < NN; k++)
        amalgCoordView[k] = coordView[k * NDim];
    }

    // On the finest level, we must map userCpts (which corresponds to
    // pressure cpts and pressure mid-points) to the velocity variables
    //
    // NOTE: on coarser levels the lower numbered velocity dofs correspond
    // to points that are co-located with pressures and the two numberings
    // are identical so no translation is needed.
    if (fineLevelID == 0) {
      ArrayRCP<LO> p2vMap = fineLevel.Get<ArrayRCP<LO> >("p2vMap", NoFactory::get());

      for (int k = 0; k < userCpts.size(); k++)
        userCpts[k] = p2vMap[userCpts[k]] / NDim;
    }

    GetOStream(Runtime1) << "Amalgamated velocity C-points: " << userCpts.size() << " " << userCpts << std::endl;

    // Now determine velocity CPOINTs for amalgamated system
    RCP<MyCptList> amalgCpts = rcp(new MyCptList(NN));
    std::vector<char> amalgStatus(NN, UNASSIGNED);

    FindDist4Cpts(*amalgA, *amalgCoords, userCpts, amalgStatus, *amalgCpts, fineLevelID);

    if (pL.get<bool>("phase2")) {
      // Beef up any limited pattern
      PhaseTwoPattern(*amalgA, *amalgCoords, amalgStatus, *amalgCpts);
    }

    // Unamalgamate data
    Array<LO>& Cptlist               = myCpts->getCList();
    Array<LO>& amalgCptlist          = amalgCpts->getCList();
    std::vector<short>& numCpts      = myCpts->getNumCpts();
    std::vector<short>& amalgNumCpts = amalgCpts->getNumCpts();

    int p = amalgCptlist.size();

    Cptlist.resize(p * NDim);
    for (int k = 0; k < p; k++) {
      Cptlist[k * NDim] = amalgCptlist[k] * NDim;
      for (int j = 1; j < NDim; j++)
        Cptlist[k * NDim + j] = Cptlist[k * NDim] + j;
    }

    for (Xpetra::global_size_t i = 0; i < NN; i++) {
      for (int j = 0; j < NDim; j++) {
        status[N - 1 - (i * NDim + j)] = amalgStatus[NN - 1 - (i)];
        numCpts[i * NDim + j]          = amalgNumCpts[i];
        for (int k = 0; k < amalgNumCpts[i]; k++)
          (*myCpts)(i * NDim + j)[k] = (*amalgCpts)(i)[k] * NDim + j;
      }
    }
  }

  const bool doStatusOutput = pL.get<bool>("dump status");
  if (doStatusOutput) {
    const Array<LO>& Cptlist    = myCpts->getCList();
    std::vector<short>& numCpts = myCpts->getNumCpts();

    std::string depPrefix = std::string("dep0-l") + toString(fineLevel.GetLevelID()) + (pressureMode ? "-p-" : "-v-");

    std::vector<char> depStatus(N);
    // Graph is unamalgamted, so we need to skip some CPOINTs as they are
    // essentially duplicated for different velocities
    for (int k = 0; k < Cptlist.size(); k += NDim) {
      for (Xpetra::global_size_t i = 0; i < N; i++) {
        bool isPresent = false;
        for (int j = 0; j < numCpts[i]; j++)
          if ((*myCpts)(i)[j] == Cptlist[k])
            isPresent = true;
        depStatus[i] = (isPresent ? FPOINT : UNASSIGNED);
      }
      depStatus[Cptlist[k]] = CPOINT;

      DumpStatus(depPrefix + toString(k), depStatus, NDim, false);
    }
  }

  RCP<Matrix> P;
  // FIXME :hardwired hack, pressure gids must not overlap with velocity gids
  if (pressureMode)
    CptDepends2Pattern(*A, *myCpts, P, 999999);
  else
    CptDepends2Pattern(*AForPat, *myCpts, P, 0);

#if 0
    if (pressureMode) {
      IO::Write("Ap_l"      + MueLu::toString(fineLevel.GetLevelID())   + ".mm", *A);
      IO::Write("Pp_tent_l" + MueLu::toString(coarseLevel.GetLevelID()) + ".mm", *P);
    } else {
      IO::Write("Av_l"      + MueLu::toString(fineLevel.GetLevelID())   + ".mm", *A);
      IO::Write("Pv_tent_l" + MueLu::toString(coarseLevel.GetLevelID()) + ".mm", *P);
    }
#endif

  // Construct coarse map
  RCP<const Map> coarseMap = P->getDomainMap();

  // Construct coarse nullspace
  RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, 1);
  coarseNullspace->putScalar(STS::one());

  // Construct coarse coordinates
  const Array<LO>& Cptlist      = myCpts->getCList();
  RCP<MultiVector> coarseCoords = MultiVectorFactory::Build(coarseMap, NDim);
  for (int k = 0; k < NDim; k++) {
    ArrayRCP<const SC> coords1D = coords->getData(k);
    ArrayRCP<SC> coarseCoords1D = coarseCoords->getDataNonConst(k);

    for (int i = 0; i < coarseCoords1D.size(); i++)
      coarseCoords1D[i] = coords1D[Cptlist[i]];
  }

  // Level Set
  Set(coarseLevel, "P", P);
  Set(fineLevel, "CoarseMap", coarseMap);
  if (pressureMode) {
    Set(coarseLevel, "CoordinatesPressure", coarseCoords);

  } else {
    Set(coarseLevel, "CoordinatesVelocity", coarseCoords);
    // FIXME: why does coarse pattern matrix look like?
    RCP<Matrix> AP  = Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*AForPat, false, *P, false, GetOStream(Statistics2), true, true);
    RCP<Matrix> RAP = Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*P, true, *AP, false, GetOStream(Statistics2), true, true);
    Set(coarseLevel, "AForPat", RAP);
  }
  Set(coarseLevel, "Nullspace", coarseNullspace);

  // Compute data for velocity
  if (pressureMode) {
    Array<LO> velCptlist = Cptlist;
    FindMidPoints(*A, *coords, velCptlist, *myCpts);
    coarseLevel.Set<Array<LO> >("PresCptsAndMids", velCptlist, NoFactory::get());
  }
}

template <class T>
void PrintVector(const std::vector<T>& v, const std::string& name, int n = -1) {
  std::cout << "======================" << std::endl;
  if (!name.empty())
    std::cout << "=== " << name << " ===" << std::endl;
  if (n == -1)
    n = v.size();
  for (int i = 0; i < n; i++)
    std::cout << i << ": " << v[i] << std::endl;
  std::cout << "======================" << std::endl;
}

// distance2 returns _squared_ distance
// i.e. the final sqrt calculation is not done
template <class SC>
SC distance2(const ArrayRCP<ArrayRCP<const SC> >& coords1D, int i, int j) {
  const int NDim = coords1D.size();

  SC d = Teuchos::ScalarTraits<SC>::zero();
  for (int k = 0; k < NDim; k++) {
    SC dtmp = coords1D[k][j] - coords1D[k][i];
    d += dtmp * dtmp;
  }

  return d;
}

template <int N = 10000>
std::string i2s(int i) {
  std::ostringstream os;
  if (i >= N) {
    os << i;
  } else {
    if (i < 10) os << "0";
    if (i < 100) os << "0";
    if (i < 1000) os << "0";
    if (i < 10000) os << "0";
    if (i < 100000) os << "0";
    if (i < 1000000) os << "0";
    if (i < 10000000) os << "0";
    if (i < 100000000) os << "0";
    if (i < 1000000000) os << "0";
  }
  os << i;
  return os.str();
}

// Initial fill Cptlist with a set of distance 4 points (during phase one).
// Additional Cpts are then determined looking for large gaps between the
// phase one Cpts. Candidate additional Cpts corresponds to phase one FPOINTs
// that have only 1 or 2 Cpts within a graph distance of 3 and are generally
// far (via graph or coordinate distances) from existing Cpts. We also define
// a sparsity pattern. An initial pattern is computed which basically
// includes all FPOINTs within a distance 3 from a Cpt. Additional entries
// are added to the initial sparsity pattern via PhaseTwoPattern(). These
// points correspond to Fpoints that only interpolate from 2 or less Cpts,
// are also far from existing Cpoints, and where the orientation of the
// interpolation Cpts is skewed to one side of the Fpoint (see
// PhaseTwoPattern for more details).
//
// NOTE: inefficiencies
// The main inefficiency is the sorting of the CandidateList. This is
// expensive and it is not clear how important it really is. My guess is that
// we could actually leave the CandidateList unsorted and things would work
// almost as well. We could also do some kind of hybrid where we sort the
// first 100 candidates. Use this data to create something like 10 bins and
// for any further candidates just put them in the right bin. When we need to
// choose a new Cpt from the CandidateList, just pick any one from the lowest
// populated bin. There are also potential inefficiences with respect to
// malloc(). I doubt that these are a big problem, but one could allocate
// some workspaces ahead of time to avoid the constant malloc/free cycle in
// CompDistances().
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindDist4Cpts(const Matrix& A, const MultiVector& coords, const Array<LO>& userCpts, std::vector<char>& status, MyCptList& myCpts, int levelID) const {
  int NDim       = coords.getNumVectors();
  size_t numRows = A.getLocalNumRows();

  ArrayRCP<const size_t> ia;
  ArrayRCP<const LO> ja;
  CreateCrsPointers(A, ia, ja);

  ArrayRCP<ArrayRCP<const SC> > coords1D(NDim);
  for (int k = 0; k < NDim; k++)
    coords1D[k] = coords.getData(k);

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero();

  // Initialize coordDist to some large value.
  // coordDist is an attempt to measure an average distance from a given
  // point to all the CPOINTs that it depends on. The averages are harmonic, so
  // basically the initial large coordDist will be averaged away with the 1st
  // harmonic average. The initial big value computed here is
  //     (Max(x)-Min(x))^2 + (Max(y)-Min(y))^2 + (Max(z)-Min(z))^2
  SC big = zero;
  for (int i = 0; i < NDim; i++) {
    SC dmin = *std::min_element(coords1D[i].begin(), coords1D[i].end());
    SC dmax = *std::max_element(coords1D[i].begin(), coords1D[i].end());

    big += ((dmax - dmin) * (dmax - dmin));
  }
  // MATCH_MATLAB
  std::vector<SC> coordDist(numRows, 10000 * big);

  const ParameterList& pL   = GetParameterList();
  const bool doStatusOutput = pL.get<bool>("dump status");
  const bool pressureMode   = (pL.get<std::string>("mode") == "pressure");

  // Set all Dirichlet points as Fpoints FIXME: why FPOINTs?
  // However, if a Dirichlet point is in userCpts, it will be added to the
  // Cpt list later
  for (size_t i = 0; i < numRows; i++)
    if (ia[i + 1] - ia[i] == 1)
      status[i] = FPOINT;

  // userCpts have already been fixed to be CPOINTs so we want to first mark
  // them appropriately, and put them first in the CPOINT list, but still go
  // through loops below to update distances and FPOINTs. Initialization is
  // done here so that these points do not end up with more than 1 nnz in the
  // associated sparsity pattern row.
  TEUCHOS_TEST_FOR_EXCEPTION(myCpts.getCList().size(), Exceptions::RuntimeError,
                             "myCpts in FindDist4Points must be uninitialized");
  Array<LO>& Cptlist = myCpts.getCList();
  for (int i = 0; i < userCpts.size(); i++) {
    status[userCpts[i]] = CPOINT_U;
    Cptlist.push_back(userCpts[i]);
  }

  std::string st = std::string("status-l") + toString(levelID) + (pressureMode ? "-p-" : "-v-");
  int dumpCount  = 0;
  if (doStatusOutput) {
    DumpCoords("coord-l" + toString(levelID) + (pressureMode ? "-p" : "-v"), coords);
    DumpStatus(st + i2s(dumpCount++) + "-A", status, NDim);
  }

  std::vector<short>& numCpts = myCpts.getNumCpts();

  // Determine CPOINTs
  int userCcount       = 0;
  size_t numCandidates = 0;
  std::vector<char> distIncrement(numRows, 0);
  std::vector<double> cumGraphDist(numRows, 0.);
  std::vector<LO> candidateList(numRows, 0);
  size_t i = 0;
  while (i < numRows) {
    LO newCpt = -1;

    // Check userCpts list
    //
    // These were essentially already determined to be CPOINTs by some other
    // functions. We want to then add them one-by-one, updating distances,
    // and look to see if further CPOINTs can be added once we have finished
    // all of the userCpts
    if (userCcount < userCpts.size())
      newCpt = userCpts[userCcount++];

    // Check for possible CPOINT on candidate list
    // FIXME: Could CANDIDATE list contain non-CANDIDATE statuses?
    while ((newCpt == -1) && (numCandidates > 0)) {
      if (status[candidateList[numCandidates - 1]] <= CANDIDATE) {
        newCpt         = candidateList[numCandidates - 1];
        status[newCpt] = CPOINT;

        if (doStatusOutput)
          DumpStatus(st + i2s(dumpCount++) + "-B", status, NDim);
      }
      numCandidates--;
      // FIXME: Why is there no i++ here?
    }

    // If no new CPOINT identified in candidate list, check the unassigned list
    while ((newCpt == -1) && (i < numRows)) {
      if (status[i] == UNASSIGNED) {
        newCpt         = i;
        status[newCpt] = CPOINT;

        if (doStatusOutput)
          DumpStatus(st + i2s(dumpCount++) + "-C", status, NDim);
      }
      i++;
    }

    // Update distances and the status of neighbors neighbors to reflect a
    // newly found CPOINT
    if (newCpt != -1) {
      std::vector<LO> dist1, dist2, dist3, dist4;
      // FIXME: Should CompDistances automatically exclude other CPOINTs?
      CompDistances(A, newCpt, 4, dist1, dist2, dist3, dist4);

      // Make sure that the only CPOINT in dist3 is newCpt. All others should be excluded.
      int numDist3 = 0;
      for (size_t k = 0; k < dist3.size(); k++) {
        LO j = dist3[k];
        if (status[j] < CPOINT)
          dist3[numDist3++] = j;
      }
      dist3.resize(numDist3);
      dist3.push_back(newCpt);

      // UNASSIGNED or CANDIDATE distance 3 and closer neighbors are put into FPOINT list
      // FIXME: why not put TWOTIMER there too?
      bool dumpStatus = false;
      for (size_t k = 0; k < dist3.size(); k++) {
        LO j = dist3[k];
        if (status[j] == UNASSIGNED || status[j] == CANDIDATE) {
          status[j]  = FPOINT;
          dumpStatus = true;
        }
      }
      if (dumpStatus && doStatusOutput)
        DumpStatus(st + i2s(dumpCount++) + "-D", status, NDim);

      // Update myCpts() to reflect dependence of neighbors on newCpt
      for (size_t k = 0; k < dist3.size(); k++) {
        LO j = dist3[k];

        TEUCHOS_TEST_FOR_EXCEPTION(numCpts[j] >= myCpts.getNnzPerRow(), Exceptions::RuntimeError,
                                   "Increase max number of C points per row");
        myCpts(j)[numCpts[j]++] = newCpt;
      }

      // Update cumGraphDist
      // NOTE: order matters as dist2 is contained within dist3, etc.
      // FIXME: Do dist2 and dist1 contain CPOINTs?
      for (size_t k = 0; k < dist3.size(); k++) distIncrement[dist3[k]] = 3;
      for (size_t k = 0; k < dist2.size(); k++) distIncrement[dist2[k]] = 2;
      for (size_t k = 0; k < dist1.size(); k++) distIncrement[dist1[k]] = 1;
      distIncrement[newCpt] = 0;

      for (size_t k = 0; k < dist3.size(); k++) {
        LO j = dist3[k];
        // MATCH_MATLAB: (numCpts[j]-1) is to match Matlab, where numCpts is updated after distance calculation
        cumGraphDist[j] = (cumGraphDist[j] * (numCpts[j] - 1) + distIncrement[j]) / numCpts[j];
      }
      cumGraphDist[newCpt] = 0;

      // Compute coordinate distance to CPOINT
      //
      // Distance of CANDIDATEs to CPOINTs will be used to determine the next
      // chosen CPOINT from the candidate list. Distances are also used to
      // decide where new CPOINTs should be added.
      for (size_t k = 0; k < dist4.size(); k++) {
        LO j = dist4[k];

        SC distance = distance2(coords1D, newCpt, j);

        // Harmonic average new distance with old distance
        // There should really be a '2' in front of this expression. This
        // is actually a bug in the code. However, if I put a '2', I don't
        // get the perfect coarsening for a uniform mesh ... so I'm leaving
        // if for now without the 2.
        // MATCH_MATLAB
        coordDist[j] = 2.0 * (coordDist[j] * distance) / (coordDist[j] + distance);
#if 0
          SC kkk = 10.;
          if (coordDist[j] > distance)
            coordDist[j] = (kkk*coordDist[j]*distance)/(coordDist[j]*(kkk-1)+ distance);
          coordDist[j] = (kkk*coordDist[j]*distance)/(coordDist[j]        + distance*(kkk-1));
#endif
      }

      // Mark all unassigned dist4 points as CANDIDATE and compress
      // dist4 so that it only contains entries for the candidate list.
      size_t numNewCandidates = 0;
      dumpStatus              = false;
      for (size_t k = 0; k < dist4.size(); k++) {
        LO j = dist4[k];

        // NOTE: numNewCandidates is always <= k, so we don't overwrite the
        // dist4 before reading
        if (status[j] == CANDIDATE) {
          // Mark as already being assigned again to candidate list so that
          // entry in old 'sorted' candidate list can be removed and a new
          // 'unsorted' entry can be created. This new entry will later be
          // sorted reflecting the new coordinate distance.
          status[j]                 = TWOTIMER;
          dist4[numNewCandidates++] = j;
          dumpStatus                = true;

        } else if (status[j] == UNASSIGNED) {
          status[j]                 = CANDIDATE;
          dist4[numNewCandidates++] = j;
          dumpStatus                = true;
        }
      }
      dist4.resize(numNewCandidates);
      if (dumpStatus && doStatusOutput)
        DumpStatus(st + i2s(dumpCount++) + "-E", status, pressureMode);

      // Now remove all TWOTIMERs from the old candidate list
      size_t numOldCandidates = 0;
      dumpStatus              = false;
      for (size_t k = 0; k < numCandidates; k++) {
        LO j = candidateList[k];

        if (status[j] == CANDIDATE) {
          candidateList[numOldCandidates++] = j;
        }
        if (status[j] == TWOTIMER) {
          status[j]  = CANDIDATE;
          dumpStatus = true;
        }
      }
      if (dumpStatus && doStatusOutput)
        DumpStatus(st + i2s(dumpCount++) + "-F", status, NDim);

      // Sort the candidates based on distances (breaking ties via degrees,
      // encouraging points near boundary). First, we order new candidates
      // and then we merge together two sorted lists.
      //
      // NOTE: to match matlab (and break ties), I added the  1.e-10 term
      std::vector<double> ddtemp(numNewCandidates);
      for (size_t k = 0; k < numNewCandidates; k++) {
        LO j = dist4[k];
        // MATCH_MATLAB
#ifdef optimal
        // This one is better, but we are trying to replicate Matlab now
        ddtemp[k] = -coordDist[j] - .01 * (ia[j + 1] - ia[j]) + 1e-10 * (j + 1);
#else
        ddtemp[k] = +coordDist[j] - 0.0 * (ia[j + 1] - ia[j]) + 1e-3 * (j + 1);
#endif
      }
      Muelu_az_dsort2(ddtemp, dist4);

      MergeSort(candidateList, numOldCandidates, dist4, coordDist, ia);
      numCandidates = numOldCandidates + numNewCandidates;
    }
  }

  // Add additional CPOINTs based on some score which includes the number of
  // CPOINTs that an FPOINT depends on as well as its distance (both graph
  // and coordinate) to nearby CPOINTs.
  const double graphWeight  = 0.8;
  const double orientWeight = 0.5;
  for (int numCDepends = 1; numCDepends <= 2; numCDepends++) {
    numCandidates = 0;

    std::vector<int> candidates;
    for (i = 0; i < numRows; i++)
      if (status[i] < CPOINT && numCpts[i] == numCDepends) {
        candidates.push_back(i);
        numCandidates++;
      }

    if (numCandidates != 0) {
      // Sort FPOINTs based on distance to CPOINTs and orientation
      double maxGraphDist = -1e20;
      double maxCoordDist = -1e20;
      for (size_t p = 0; p < numCandidates; p++) {
        LO j = candidates[p];

        maxGraphDist = std::max(maxGraphDist, cumGraphDist[j]);
        maxCoordDist = std::max(maxCoordDist, coordDist[j]);
      }

      std::vector<double> score(numCandidates);
      std::vector<double> orientation(numCandidates);
      for (size_t p = 0; p < numCandidates; p++) {
        LO j = candidates[p];

        double graphScore = cumGraphDist[j] / maxGraphDist;
        double coordScore = coordDist[j] / maxCoordDist;
        // MATCH_MATLAB
        score[p] = -(graphWeight * graphScore + (1 - graphWeight) * coordScore + 1e-6 * (j + 1));

        if (numCDepends == 2) {
          // Orientation of -1 means that we have two CPOINTs on opposite
          // sides of the FPOINT. Orientation of 0 implies that things are
          // orthogonal. orientation of 1 implies that Cpoints are on the
          // same side
          SC norm = zero, vec1[3], vec2[3];
          for (int k = 0; k < NDim; k++) {
            vec1[k] = coords1D[k][j] - coords1D[k][myCpts(j)[0]];
            norm += vec1[k] * vec1[k];
          }
          norm = sqrt(norm);
          for (int k = 0; k < NDim; k++)
            vec1[k] /= norm;

          norm = zero;
          for (int k = 0; k < NDim; k++) {
            vec2[k] = coords1D[k][j] - coords1D[k][myCpts(j)[1]];
            norm += vec2[k] * vec2[k];
          }
          norm = sqrt(norm);
          for (int k = 0; k < NDim; k++)
            vec2[k] /= norm;

          orientation[p] = zero;
          for (int k = 0; k < NDim; k++)
            orientation[p] += vec1[k] * vec2[k];

          score[p] = -(orientWeight * orientation[p] - (1 - orientWeight) * score[p]);

        } else {
          orientation[p] = 1.0;
        }
      }

      std::vector<LO> index(numCandidates);
      for (size_t p = 0; p < numCandidates; p++)
        index[p] = p;
      Muelu_az_dsort2(score, index);

      for (size_t p = 0; p < numCandidates; p++) {
        int newCpt = candidates[index[p]];

        if (numCpts[newCpt] == numCDepends &&
            cumGraphDist[newCpt] >= 2.6 &&
            orientation[index[p]] > -0.2) {
          status[newCpt]    = CPOINT;
          numCpts[newCpt]   = 1;
          myCpts(newCpt)[0] = newCpt;

          if (doStatusOutput)
            DumpStatus(st + i2s(dumpCount++) + "-G", status, NDim);

          std::vector<LO> dist1, dist2, dist3, dist4;
          CompDistances(A, newCpt, 3, dist1, dist2, dist3, dist4);

          // Make sure that there are no CPOINTs in dist1, dist2, dist3.
          int numDist1 = 0;
          for (size_t k = 0; k < dist1.size(); k++) {
            LO j = dist1[k];
            if (status[j] < CPOINT)
              dist1[numDist1++] = j;
          }
          dist1.resize(numDist1);
          int numDist2 = 0;
          for (size_t k = 0; k < dist2.size(); k++) {
            LO j = dist2[k];
            if (status[j] < CPOINT)
              dist2[numDist2++] = j;
          }
          dist2.resize(numDist2);
          int numDist3 = 0;
          for (size_t k = 0; k < dist3.size(); k++) {
            LO j = dist3[k];
            if (status[j] < CPOINT)
              dist3[numDist3++] = j;
          }
          dist3.resize(numDist3);

          // Update cumGraphDist
          // NOTE: order matters as dist2 is contained within dist3, etc.
          for (size_t k = 0; k < dist3.size(); k++) distIncrement[dist3[k]] = 3;
          for (size_t k = 0; k < dist2.size(); k++) distIncrement[dist2[k]] = 2;
          for (size_t k = 0; k < dist1.size(); k++) distIncrement[dist1[k]] = 1;
          distIncrement[newCpt] = 0;

          // Update myCpts() to reflect dependence of neighbors on newCpt
          for (size_t k = 0; k < dist3.size(); k++) {
            LO j = dist3[k];

            TEUCHOS_TEST_FOR_EXCEPTION(numCpts[j] >= myCpts.getNnzPerRow(), Exceptions::RuntimeError,
                                       "Increase max number of C points per row");
            myCpts(j)[numCpts[j]++] = newCpt;
          }

          for (size_t k = 0; k < dist3.size(); k++) {
            LO j = dist3[k];
            // (numCpts[j]-1) is to match Matlab, where numCpts is updated after distance calculation
            cumGraphDist[j] = (cumGraphDist[j] * (numCpts[j] - 1) + distIncrement[j]) / numCpts[j];
          }
          cumGraphDist[newCpt] = 0;
        }
      }
    }
  }

  // Build up the CPOINT list
  for (i = 0; i < numRows; i++)
    if (status[i] == CPOINT)
      Cptlist.push_back(i);
    else if (status[i] == CPOINT_U)
      status[i] = CPOINT;
}

// Look at pattern rows which have only one or two nonzeros and consider
// adding additional nonzero entries. New nonzero possibilities for row k are
// obtained by looking at k's neighbors (as determined by the matrix) to see
// what CPOINTs these neighbors interpolate from.  The final determination is
// based on a composite score that considers the distance between k and the
// possible new CPOINT as well as the orientation of the new possible CPOINT
// with respect to k's current CPOINTs. Generally, points which are on the
// opposite side of k' current CPOINTs are favored.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    PhaseTwoPattern(const Matrix& A, const MultiVector& coords, const std::vector<char>& status, MyCptList& myCpts) const {
  GetOStream(Runtime0) << "Starting phase 2" << std::endl;

  int NDim       = coords.getNumVectors();
  size_t numRows = A.getLocalNumRows();

  ArrayRCP<const size_t> ia;
  ArrayRCP<const LO> ja;
  CreateCrsPointers(A, ia, ja);

  ArrayRCP<ArrayRCP<const SC> > coords1D(NDim);
  for (int k = 0; k < NDim; k++)
    coords1D[k] = coords.getData(k);

  std::vector<short>& numCpts = myCpts.getNumCpts();

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero();

  size_t N = myCpts.getCList().size();
  std::vector<int> nearbyCs(N);
  std::vector<double> score(N);
  std::vector<double> dists(N);

  std::vector<char> scratch(numRows, 'n');
  std::vector<LO> candidates(numRows);

  for (int numCDepends = 1; numCDepends <= 2; numCDepends++) {
    int numCandidates = 0;
    for (size_t i = 0; i < numRows; i++)
      if (numCpts[i] == numCDepends && status[i] < CPOINT)
        candidates[numCandidates++] = i;

    for (int p = 0; p < numCandidates; p++) {
      // Mark already existing CPOINT dependencies
      LO i     = candidates[p];
      LO* cpts = myCpts(i);
      for (int k = 0; k < numCpts[i]; k++)
        scratch[cpts[k]] = 'y';

      // Make a list of my neighbors' CPOINT dependencies, excluding all
      // already existing CPOINT dependencies for candidates[p]
      const LO* neighs = &ja[ia[i]];
      int numNeighbors = ia[i + 1] - ia[i];
      int numNearbyCs  = 0;
      for (int k = 0; k < numNeighbors; k++) {
        LO curNeigh       = neighs[k];
        const LO* neighCs = myCpts(curNeigh);

        for (int j = 0; j < numCpts[curNeigh]; j++) {
          LO neighNeighC = neighCs[j];

          if (scratch[neighNeighC] != 'y') {
            scratch[neighNeighC]    = 'y';
            nearbyCs[numNearbyCs++] = neighNeighC;
          }
        }
      }

      // Reset scratch
      for (int k = 0; k < numCpts[i]; k++)
        scratch[cpts[k]] = 'n';
      for (int k = 0; k < numNearbyCs; k++)
        scratch[nearbyCs[k]] = 'n';

      // MATCH_MATLB
      std::sort(nearbyCs.begin(), nearbyCs.begin() + numNearbyCs);

      if (numNearbyCs != 0) {
        SC norm = zero, vec1[3] = {0.0, 0.0, 0.0}, vec2[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < NDim; k++) {
          vec1[k] = coords1D[k][i] - coords1D[k][cpts[0]];
          norm += vec1[k] * vec1[k];
        }
        norm = sqrt(norm);
        for (int k = 0; k < NDim; k++)
          vec1[k] /= norm;

        if (numCDepends == 2) {
          norm = zero;
          for (int k = 0; k < NDim; k++) {
            vec2[k] = coords1D[k][i] - coords1D[k][cpts[1]];
            norm += vec2[k] * vec2[k];
          }
          norm = sqrt(norm);
          for (int k = 0; k < NDim; k++)
            vec2[k] /= norm;

        } else {
          for (int k = 0; k < NDim; k++)
            vec2[k] = vec1[k];
        }

        for (int k = 0; k < numNearbyCs; k++) score[k] = 0;
        for (int k = 0; k < numNearbyCs; k++) dists[k] = 0;

        for (int j = 0; j < numNearbyCs; j++) {
          SC newVec[3];

          norm = 0;
          for (int k = 0; k < NDim; k++) {
            newVec[k] = coords1D[k][nearbyCs[j]] - coords1D[k][i];
            norm += newVec[k] * newVec[k];
          }
          norm = sqrt(norm);
          for (int k = 0; k < NDim; k++)
            newVec[k] /= norm;

          score[j] = 0;
          for (int k = 0; k < NDim; k++)
            score[j] += newVec[k] * (vec1[k] + vec2[k]);
          // Why??
          score[j] /= 2;

          dists[j] = norm;
        }

        // Normalize distances
        double maxDist = 0.;
        for (int j = 0; j < numNearbyCs; j++)
          if (maxDist < dists[j])
            maxDist = dists[j];

        for (int j = 0; j < numNearbyCs; j++)
          dists[j] /= maxDist;

        const double distWeight = 0.3;
        double maxComposite     = -10000;
        int maxIndex            = -1;
        for (int j = 0; j < numNearbyCs; j++) {
          // The formula is
          //     if (score[j] - distWeight*dists[j] > maxComposite)
          // MATCH_MATLAB
          double composite = score[j] - distWeight * dists[j] + 1.0e-7 * (nearbyCs[j] - 1);
          if (maxComposite < composite) {
            maxComposite = composite;
            maxIndex     = j;
          }
        }

        if (score[maxIndex] - 0.2 * numCDepends > -0.3) {
          TEUCHOS_TEST_FOR_EXCEPTION(numCpts[i] >= myCpts.getNnzPerRow(),
                                     Exceptions::RuntimeError, "Increase max number of C points per row");
          myCpts(i)[numCpts[i]++] = nearbyCs[maxIndex];
        }
      }
    }
  }
}

// Compute mid-points associated with adjacent Cpts.
//
// Basically, the location of each Fpoint is compared with the average
// location of each Cpoint in Pat(Fpoint,:). If the location of the Fpoint is
// close to this average it is considered as a possible mid-point. In
// addition, however, we look to see if a possible mid-point is "close" or
// not to an already computed mid-point. If it is NOT too close, then this
// possible mid-point is declared to be an actual mid-point.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindMidPoints(const Matrix& A, const MultiVector& coords, Array<LO>& Cptlist, const MyCptList& myCpts) const {
  int NDim       = coords.getNumVectors();
  size_t numRows = A.getLocalNumRows();

  const ParameterList& pL = GetParameterList();
  double tau_2            = pL.get<double>("tau_2");
  // In calculations, we use tau_2^2
  tau_2 = tau_2 * tau_2;

  const std::vector<short>& numCpts = myCpts.getNumCpts();

  ArrayRCP<const size_t> ia;
  ArrayRCP<const LO> ja;
  CreateCrsPointers(A, ia, ja);

  ArrayRCP<ArrayRCP<const SC> > coords1D(NDim);
  for (int k = 0; k < NDim; k++)
    coords1D[k] = coords.getData(k);

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero();

  // Calculate number of nonzeros per row, make it negative, and then sort.
  // The idea is that when assigning midpoints, we want to start by looking
  // at points which have many coarse point dependencies
  std::vector<int> nnzPerRow(numRows);
  std::vector<int> index(numRows);
  for (size_t i = 0; i < numRows; i++) {
    nnzPerRow[i] = (numCpts[i] ? numCpts[i] : 1);
    nnzPerRow[i] = -100000 * numCpts[i] + i;
    index[i]     = i;
  }

  // Sort only for the purposes of filling 'index', which determines the
  // order that we search for possible midpoints
  Muelu_az_sort<SC>(&nnzPerRow[0], numRows, &index[0], NULL);

  // Reset so that we have unsorted version of nnzPerRow and also mark points
  // which cannot be mid points
  std::vector<char> lookedAt(numRows, 'n');
  for (size_t i = 0; i < numRows; i++) {
    nnzPerRow[i] = (numCpts[i] ? numCpts[i] : 1);
    if (nnzPerRow[i] == 1)
      lookedAt[i] = 'y';
  }
  for (int i = 0; i < Cptlist.size(); i++)
    lookedAt[Cptlist[i]] = 'y';

  // Compute some target midpoints based on taking averages associated with
  // the sparsity pattern and coarse grid point locations.
  ArrayRCP<ArrayRCP<SC> > targetMidCoords1D(NDim);
  for (int k = 0; k < NDim; k++) {
    ArrayRCP<SC>& target1D = targetMidCoords1D[k];

    target1D.resize(numRows);
    for (size_t i = 0; i < numRows; i++) {
      target1D[i] = zero;

      for (int j = 0; j < numCpts[i]; j++)
        target1D[i] += coords1D[k][myCpts(i)[j]];

      target1D[i] /= nnzPerRow[i];
    }
  }

  std::vector<char> isMidPoint(numRows, 'n');
  std::vector<char> inNearbyCs(numRows, 'n');
  std::vector<char> inNeighs(numRows, 'n');
  std::vector<int> neighs(numRows);
  std::vector<int> sameCGroup(50);

  int numMidPoints = 0;
  for (size_t i = 0; i < numRows; i++) {
    int curF = index[i];

    if (lookedAt[curF] == 'y')
      continue;
    lookedAt[curF] = 'y';

    const LO* curFCs = myCpts(curF);

    for (int j = 0; j < numCpts[curF]; j++)
      inNearbyCs[curFCs[j]] = 'y';

    // Find all FPOINTs with the same curFCs (perhaps
    // containing additional curFCs) as curF and
    // put them in sameCGroup
    int numNeigh       = 0;
    neighs[numNeigh++] = curF;
    inNeighs[curF]     = 'y';

    int nextLayerStart = 0;
    int nextLayerEnd   = 0;
    int numSameGrp     = 0;

    int flag = 1;
    while (flag == 1) {
      flag = 0;

      for (int k = nextLayerStart; k <= nextLayerEnd; k++) {
        LO curNeigh       = neighs[k];
        const LO* neighCs = myCpts(curNeigh);

        // Check if subset of this neighbor's CPOINT dependencies include all
        // the CPOINT dependencies of curF
        int sum = 0;
        for (int j = 0; j < numCpts[curNeigh]; j++)
          if (inNearbyCs[neighCs[j]] == 'y')
            sum++;

        if (sum == nnzPerRow[curF]) {
          lookedAt[curNeigh] = 'y';

          // Make sure we have enough space
          if (Teuchos::as<int>(sameCGroup.size()) <= numSameGrp)
            sameCGroup.resize(2 * numSameGrp);

          sameCGroup[numSameGrp++] = curNeigh;
          flag                     = 1;
        }

        // Add neighbors of curNeigh that haven't already been
        // add to the neighbor list while processing curF
        for (size_t j = ia[curNeigh]; j < ia[curNeigh + 1]; j++)
          if (inNeighs[ja[j]] == 'n') {
            neighs[numNeigh++] = ja[j];
            inNeighs[ja[j]]    = 'y';
          }
      }

      nextLayerStart = nextLayerEnd + 1;
      nextLayerEnd   = numNeigh - 1;
    }

    // Reset status arrays
    for (int j = 0; j < numNeigh; j++)
      inNeighs[neighs[j]] = 'n';
    for (int j = 0; j < numCpts[curF]; j++)
      inNearbyCs[curFCs[j]] = 'n';

    // At this point we have now constructed a group of possible mid points
    // all with the same Cpt dependencies. Now, we need to find the one in
    // this group which is closest to the target midpoint coordinates.
    double smallest   = 1e30;
    int smallestIndex = -1;
    for (int j = 0; j < numSameGrp; j++) {
      // MATCH_MATLAB
      double dist = 1e-8 * (sameCGroup[j] + 1);

      for (int k = 0; k < NDim; k++) {
        double dtemp = coords1D[k][sameCGroup[j]] - targetMidCoords1D[k][curF];
        dist += dtemp * dtemp;
      }
      if (dist < smallest) {
        smallest      = dist;
        smallestIndex = sameCGroup[j];
      }
    }

    // So now smallestIndex is the best midpoint candidate within sameCGroup.
    // We now need to check if smallestIndex is really close to an already
    // existing mid-point. In fact, we could have multiple copies of
    // mid-points or some very close mid-points. To see this, consider
    //
    //                    P1              P2
    //
    //
    //
    //                    P3              P4
    // where P1/P4's midpoint is the same as P2/P3's midpoint. We get rid of
    // these by checking if a new potential midpoint is close to any previous
    // midpoints.

    // Check if anybody in sameCGroup is already a midpoint. If so, check
    // each one of these midpoints to see if any of these is real close to
    // the curF.
    flag = 0;
    for (int j = 0; j < numSameGrp; j++)
      if (isMidPoint[sameCGroup[j]] == 'y')
        flag = 1;

    if (flag == 1) {
      // Get an idea of the spacing between curFCs
      double delta = 0.0;
      for (int k = 0; k < NDim; k++) {
        double dmin = coords1D[k][curFCs[0]];
        double dmax = dmin;

        for (int j = 1; j < numCpts[curF]; j++) {
          SC c = coords1D[k][curFCs[j]];
          if (c < dmin) dmin = c;
          if (c > dmax) dmax = c;
        }
        delta += ((dmax - dmin) * (dmax - dmin));
      }

      // Now find the closest point among all sameCGroup midPoints to
      // smallestIndex. If this point is not too close, we make smallestIndex
      // a new mid-point.
      double close = 1000000.;
      for (int j = 0; j < numSameGrp; j++) {
        int t = sameCGroup[j];
        if (isMidPoint[t] == 'y') {
          double current = distance2(coords1D, smallestIndex, t);
          if (current < close)
            close = current;
        }
      }

      if (close / delta > tau_2) {
        isMidPoint[smallestIndex] = 'y';
        numMidPoints++;
      }

    } else {
      isMidPoint[smallestIndex] = 'y';
      numMidPoints++;
    }
  }

  // This loop also sorts mid points
  int count = 0;
  for (size_t i = 0; i < numRows; i++)
    if (isMidPoint[i] == 'y') {
      Cptlist.push_back(i);
      count++;
    }

  TEUCHOS_TEST_FOR_EXCEPTION(count != numMidPoints, Exceptions::RuntimeError,
                             "Wrong with the number of mid points: " << count << " vs. " << numMidPoints);
}

// Convert information in Cptlist, myCpts into a sparsity pattern matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CptDepends2Pattern(const Matrix& A, const MyCptList& myCpts, RCP<Matrix>& P, LocalOrdinal offset) const {
  RCP<const Map> rowMap = A.getRowMap();
  size_t numRows        = myCpts.getLocalNumRows();

  // FIXME: how does offset play here?
  const Array<LO>& Cptlist = myCpts.getCList();
  RCP<const Map> coarseMap = MapFactory::Build(rowMap->lib(), Cptlist.size(), rowMap->getIndexBase() + offset, rowMap->getComm());

  P                   = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0));
  RCP<CrsMatrix> Pcrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  ArrayRCP<size_t> iaP;
  ArrayRCP<LO> jaP;
  ArrayRCP<SC> valP;

  const std::vector<short>& numCpts = myCpts.getNumCpts();
  size_t nnzEstimate                = std::accumulate(numCpts.begin(), numCpts.end(), 0);

  Pcrs->allocateAllValues(nnzEstimate, iaP, jaP, valP);

  ArrayView<size_t> ia = iaP();
  ArrayView<LO> ja     = jaP();
  ArrayView<SC> val    = valP();

  std::vector<GO> coarseCmap(numRows, -1);
  for (int i = 0; i < Cptlist.size(); i++)
    coarseCmap[Cptlist[i]] = i;

  SC one = Teuchos::ScalarTraits<SC>::one();

  // Build up the prolongator sparsity pattern and the initial
  // guess used for Emin (which must have row sums equal to one)

  ia[0]           = 0;
  size_t nnzCount = 0;
  for (size_t i = 0; i < numRows; i++) {
    const LO* cpts = myCpts(i);

    for (int j = 0; j < numCpts[i]; j++) {
      ja[nnzCount] = coarseCmap[cpts[j]];
      TEUCHOS_TEST_FOR_EXCEPTION(ja[nnzCount] == -1, Exceptions::RuntimeError,
                                 "Point " << cpts[j] << " is not in the global list of cpoints, but it is in the list for " << i);
      val[nnzCount] = one / ((SC)numCpts[i]);
      nnzCount++;
    }
    // NOTE: we could theoretically sort here
    // Do we need to, though?
    ia[i + 1] = nnzCount;
  }

  if (rowMap->lib() == Xpetra::UseTpetra) {
    // - Cannot resize for Epetra, as it checks for same pointers
    // - Need to resize for Tpetra, as it check ().size() == ia[numRows]
    // NOTE: these invalidate ja and val views
    jaP.resize(nnzCount);
    valP.resize(nnzCount);
  }

  Pcrs->setAllValues(iaP, jaP, valP);
  Pcrs->expertStaticFillComplete(coarseMap, A.getDomainMap());
}

// Compute all points which are within a distance 1-4 from StartPt
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CompDistances(const Matrix& A, LO start, int numDist, std::vector<LO>& dist1, std::vector<LO>& dist2, std::vector<LO>& dist3, std::vector<LO>& dist4) const {
  TEUCHOS_TEST_FOR_EXCEPTION(numDist < 1 || numDist > 4, Exceptions::InvalidArgument, "CompDistances() cannot compute " << numDist << " distances");

  size_t numRows = A.getGlobalNumRows();

  ArrayRCP<const size_t> ia;
  ArrayRCP<const LO> ja;
  CreateCrsPointers(A, ia, ja);

  std::vector<char> added(numRows, 'n');
  std::vector<LO> neighs;
  neighs.reserve(100);

  neighs.push_back(start);
  added[start] = 'y';

  for (int k = 1; k <= numDist; k++) {
    int numNeighs = neighs.size();
    for (int i = 0; i < numNeighs; i++)
      for (size_t j = ia[neighs[i]]; j < ia[neighs[i] + 1]; j++)
        if (added[ja[j]] == 'n') {
          added[ja[j]] = 'y';
          neighs.push_back(ja[j]);
        }

    if (k == 1) dist1 = neighs;
    if (k == 2) dist2 = neighs;
    if (k == 3) dist3 = neighs;
    if (k == 4) dist4 = neighs;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CreateCrsPointers(const Matrix& A, ArrayRCP<const size_t>& ia, ArrayRCP<const LO>& ja) const {
  RCP<const CrsMatrixWrap> Awrap = rcp_dynamic_cast<const CrsMatrixWrap>(rcpFromRef(A));
  TEUCHOS_TEST_FOR_EXCEPTION(Awrap.is_null(), Exceptions::RuntimeError, "A is not of CrsMatrixWrap type");

  ArrayRCP<const SC> val;
  Awrap->getCrsMatrix()->getAllValues(ia, ja, val);
}

const std::string OUTPUT_DIR = "status/";

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DumpStatus(const std::string& filename, const std::vector<char>& status, int NDim, bool isAmalgamated) const {
  const std::string dirName = OUTPUT_DIR;

  struct stat sb;
  TEUCHOS_TEST_FOR_EXCEPTION(stat(dirName.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode), Exceptions::RuntimeError,
                             "Please create a \"" << dirName << "\" directory");

  std::ofstream ofs((dirName + filename).c_str());
  size_t step = (isAmalgamated ? 1 : NDim);
  for (size_t i = 0; i < status.size(); i += step)
    ofs << status[i] << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DumpCoords(const std::string& filename, const MultiVector& coords) const {
  const std::string dirName = OUTPUT_DIR;

  struct stat sb;
  if (stat(dirName.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode))
    GetOStream(Errors) << "Please create a \"" << dirName << "\" directory" << std::endl;

  const int NDim = coords.getNumVectors();
  const int n    = coords.getLocalLength();

  ArrayRCP<ArrayRCP<const SC> > coords1D(NDim);
  for (int k = 0; k < NDim; k++)
    coords1D[k] = coords.getData(k);

  std::ofstream ofs((dirName + filename).c_str());
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < NDim; k++)
      ofs << " " << coords1D[k][i];
    ofs << std::endl;
  }
}

}  // namespace MueLu

#endif  // MUELU_Q2Q1UPFACTORY_DECL_HPP
