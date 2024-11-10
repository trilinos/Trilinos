// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LINEDETECTIONFACTORY_DEF_HPP
#define MUELU_LINEDETECTIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
//#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_LineDetectionFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("linedetection: orientation");
  SET_VALID_ENTRY("linedetection: num layers");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Generating factory for coorindates");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");

  // The factory needs the information about the number of z-layers. While this information is
  // provided by the user for the finest level, the factory itself is responsible to provide the
  // corresponding information on the coarser levels. Since a factory cannot be dependent on itself
  // we use the NoFactory class as generator class, but remove the UserData keep flag, such that
  // "NumZLayers" is part of the request/release mechanism.
  // Please note, that this prevents us from having several (independent) CoarsePFactory instances!
  // TODO: allow factory to dependent on self-generated data for TwoLevelFactories -> introduce ExpertRequest/Release in Level
  currentLevel.DeclareInput("NumZLayers", NoFactory::get(), this);
  currentLevel.RemoveKeepFlag("NumZLayers", NoFactory::get(), MueLu::UserData);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Line detection (Ray style)", currentLevel);

  LO NumZDir = 0;
  RCP<CoordinateMultiVector> fineCoords;
  ArrayRCP<coordinate_type> x, y, z;
  coordinate_type *xptr = NULL, *yptr = NULL, *zptr = NULL;

  // obtain general variables
  RCP<Matrix> A         = Get<RCP<Matrix> >(currentLevel, "A");
  LO BlkSize            = A->GetFixedBlockSize();
  RCP<const Map> rowMap = A->getRowMap();
  LO Ndofs              = rowMap->getLocalNumElements();
  LO Nnodes             = Ndofs / BlkSize;

  // collect information provided by user
  const ParameterList& pL           = GetParameterList();
  const std::string lineOrientation = pL.get<std::string>("linedetection: orientation");

  // interpret "line orientation" parameter provided by the user on the finest level
  if (currentLevel.GetLevelID() == 0) {
    if (lineOrientation == "vertical")
      Zorientation_ = VERTICAL;
    else if (lineOrientation == "horizontal")
      Zorientation_ = HORIZONTAL;
    else if (lineOrientation == "coordinates")
      Zorientation_ = GRID_SUPPLIED;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "LineDetectionFactory: The parameter 'semicoarsen: line orientation' must be either 'vertical', 'horizontal' or 'coordinates'.");
  }

  // TEUCHOS_TEST_FOR_EXCEPTION(Zorientation_!=VERTICAL, Exceptions::RuntimeError, "LineDetectionFactory: The 'horizontal' or 'coordinates' have not been tested!!!. Please remove this exception check and carefully test these modes!");

  // obtain number of z layers (variable over levels)
  // This information is user-provided on the finest level and transferred to the coarser
  // levels by the SemiCoarsenPFactor using the internal "NumZLayers" variable.
  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("NumZLayers", NoFactory::get())) {
      NumZDir = currentLevel.Get<LO>("NumZLayers", NoFactory::get());  // obtain info
      GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << " (information from Level(0))" << std::endl;
    } else {
      // check whether user provides information or it can be reconstructed from coordinates
      NumZDir = pL.get<LO>("linedetection: num layers");
      if (NumZDir == -1) {
        bool CoordsAvail = currentLevel.IsAvailable("Coordinates");

        if (CoordsAvail == true) {
          // try to reconstruct the number of layers from coordinates
          fineCoords = Get<RCP<CoordinateMultiVector> >(currentLevel, "Coordinates");
          TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3, Exceptions::RuntimeError, "Three coordinates arrays must be supplied if line detection orientation not given.");
          x    = fineCoords->getDataNonConst(0);
          y    = fineCoords->getDataNonConst(1);
          z    = fineCoords->getDataNonConst(2);
          xptr = x.getRawPtr();
          yptr = y.getRawPtr();
          zptr = z.getRawPtr();

          LO NumCoords = Ndofs / BlkSize;

          /* sort coordinates so that we can order things according to lines */
          Teuchos::ArrayRCP<LO> TOrigLoc            = Teuchos::arcp<LO>(NumCoords);
          LO* OrigLoc                               = TOrigLoc.getRawPtr();
          Teuchos::ArrayRCP<coordinate_type> Txtemp = Teuchos::arcp<coordinate_type>(NumCoords);
          coordinate_type* xtemp                    = Txtemp.getRawPtr();
          Teuchos::ArrayRCP<coordinate_type> Tytemp = Teuchos::arcp<coordinate_type>(NumCoords);
          coordinate_type* ytemp                    = Tytemp.getRawPtr();
          Teuchos::ArrayRCP<coordinate_type> Tztemp = Teuchos::arcp<coordinate_type>(NumCoords);
          coordinate_type* ztemp                    = Tztemp.getRawPtr();

          // sort coordinates in {x,y,z}vals (returned in {x,y,z}temp) so that we can order things according to lines
          // switch x and y coordinates for semi-coarsening...
          sort_coordinates(NumCoords, OrigLoc, xptr, yptr, zptr, xtemp, ytemp, ztemp, true);

          /* go through each vertical line and populate blockIndices so all   */
          /* dofs within a PDE within a vertical line correspond to one block.*/
          LO NumBlocks           = 0;
          LO NumNodesPerVertLine = 0;
          LO index               = 0;

          while (index < NumCoords) {
            coordinate_type xfirst = xtemp[index];
            coordinate_type yfirst = ytemp[index];
            LO next                = index + 1;
            while ((next != NumCoords) && (xtemp[next] == xfirst) &&
                   (ytemp[next] == yfirst))
              next++;
            if (NumBlocks == 0) {
              NumNodesPerVertLine = next - index;
            }
            // the number of vertical lines must be the same on all processors
            // TAW: Sep 14 2015: or zero as we allow "empty" processors
            // TEUCHOS_TEST_FOR_EXCEPTION(next-index != NumNodesPerVertLine,Exceptions::RuntimeError, "Error code only works for constant block size now!!!\n");
            NumBlocks++;
            index = next;
          }

          NumZDir = NumNodesPerVertLine;
          GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << " (information reconstructed from provided node coordinates)" << std::endl;
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "LineDetectionFactory: BuildP: User has to provide valid number of layers (e.g. using the 'line detection: num layers' parameter).");
        }
      } else {
        GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << " (information provided by user through 'line detection: num layers')" << std::endl;
      }
    }  // end else (user provides information or can be reconstructed) on finest level
  } else {
    // coarse level information
    // TODO get rid of NoFactory here and use SemiCoarsenPFactory as source of NumZLayers instead.
    if (currentLevel.IsAvailable("NumZLayers", NoFactory::get())) {
      NumZDir = currentLevel.Get<LO>("NumZLayers", NoFactory::get());  // obtain info
      GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << std::endl;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "LineDetectionFactory: BuildP: No NumZLayers variable found. This cannot be.");
    }
  }

  // plausibility check and further variable collection
  if (Zorientation_ == GRID_SUPPLIED) {  // On finest level, fetch user-provided coordinates if available...
    bool CoordsAvail = currentLevel.IsAvailable("Coordinates");

    if (CoordsAvail == false) {
      if (currentLevel.GetLevelID() == 0)
        throw Exceptions::RuntimeError("Coordinates must be supplied if line detection orientation not given.");
      else
        throw Exceptions::RuntimeError("Coordinates not generated by previous invocation of LineDetectionFactory's BuildP() method.");
    }
    fineCoords = Get<RCP<CoordinateMultiVector> >(currentLevel, "Coordinates");
    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3, Exceptions::RuntimeError, "Three coordinates arrays must be supplied if line detection orientation not given.");
    x    = fineCoords->getDataNonConst(0);
    y    = fineCoords->getDataNonConst(1);
    z    = fineCoords->getDataNonConst(2);
    xptr = x.getRawPtr();
    yptr = y.getRawPtr();
    zptr = z.getRawPtr();
  }

  // perform line detection
  if (NumZDir > 0) {
    LO *LayerId, *VertLineId;
    Teuchos::ArrayRCP<LO> TLayerId    = Teuchos::arcp<LO>(Nnodes);
    LayerId                           = TLayerId.getRawPtr();
    Teuchos::ArrayRCP<LO> TVertLineId = Teuchos::arcp<LO>(Nnodes);
    VertLineId                        = TVertLineId.getRawPtr();

    NumZDir = ML_compute_line_info(LayerId, VertLineId, Ndofs, BlkSize,
                                   Zorientation_, NumZDir, xptr, yptr, zptr, *(rowMap->getComm()));
    // it is NumZDir=NCLayers*NVertLines*DofsPerNode;

    // store output data on current level
    // The line detection data is used by the SemiCoarsenPFactory and the line smoothers in Ifpack/Ifpack2
    Set(currentLevel, "CoarseNumZLayers", NumZDir);
    Set(currentLevel, "LineDetection_Layers", TLayerId);
    Set(currentLevel, "LineDetection_VertLineIds", TVertLineId);
  } else {
    Teuchos::ArrayRCP<LO> TLayerId        = Teuchos::arcp<LO>(0);
    Teuchos::ArrayRCP<LO> TVertLineId     = Teuchos::arcp<LO>(0);
    Teuchos::ArrayRCP<LO> TVertLineIdSmoo = Teuchos::arcp<LO>(0);

    // store output data on current level
    // The line detection data is used by the SemiCoarsenPFactory and the line smoothers in Ifpack/Ifpack2
    Set(currentLevel, "CoarseNumZLayers", NumZDir);
    Set(currentLevel, "LineDetection_Layers", TLayerId);
    Set(currentLevel, "LineDetection_VertLineIds", TVertLineId);
  }

  // automatically switch to vertical mode on the coarser levels
  if (Zorientation_ != VERTICAL)
    Zorientation_ = VERTICAL;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_compute_line_info(LocalOrdinal LayerId[], LocalOrdinal VertLineId[], LocalOrdinal Ndof, LocalOrdinal DofsPerNode, LocalOrdinal MeshNumbering, LocalOrdinal NumNodesPerVertLine, typename Teuchos::ScalarTraits<Scalar>::coordinateType* xvals, typename Teuchos::ScalarTraits<Scalar>::coordinateType* yvals, typename Teuchos::ScalarTraits<Scalar>::coordinateType* zvals, const Teuchos::Comm<int>& /* comm */) const {
  LO Nnodes, NVertLines, MyNode;
  LO NumCoords, next;  //, subindex, subnext;
  coordinate_type xfirst, yfirst;
  coordinate_type *xtemp, *ytemp, *ztemp;
  LO* OrigLoc;
  LO i, j, count;
  LO RetVal;

  RetVal = 0;
  if ((MeshNumbering != VERTICAL) && (MeshNumbering != HORIZONTAL)) {
    if ((xvals == NULL) || (yvals == NULL) || (zvals == NULL)) RetVal = -1;
  } else {
    if (NumNodesPerVertLine == -1) RetVal = -4;
    if (((Ndof / DofsPerNode) % NumNodesPerVertLine) != 0) RetVal = -3;
  }
  if ((Ndof % DofsPerNode) != 0) RetVal = -2;

  TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -1, Exceptions::RuntimeError, "Not semicoarsening as no mesh numbering information or coordinates are given\n");
  TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -4, Exceptions::RuntimeError, "Not semicoarsening as the number of z nodes is not given.\n");
  TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -3, Exceptions::RuntimeError, "Not semicoarsening as the total number of nodes is not evenly divisible by the number of z direction nodes .\n");
  TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -2, Exceptions::RuntimeError, "Not semicoarsening as something is off with the number of degrees-of-freedom per node.\n");

  Nnodes = Ndof / DofsPerNode;
  for (MyNode = 0; MyNode < Nnodes; MyNode++) VertLineId[MyNode] = -1;
  for (MyNode = 0; MyNode < Nnodes; MyNode++) LayerId[MyNode] = -1;

  if (MeshNumbering == VERTICAL) {
    for (MyNode = 0; MyNode < Nnodes; MyNode++) {
      LayerId[MyNode]    = MyNode % NumNodesPerVertLine;
      VertLineId[MyNode] = (MyNode - LayerId[MyNode]) / NumNodesPerVertLine;
    }
  } else if (MeshNumbering == HORIZONTAL) {
    NVertLines = Nnodes / NumNodesPerVertLine;
    for (MyNode = 0; MyNode < Nnodes; MyNode++) {
      VertLineId[MyNode] = MyNode % NVertLines;
      LayerId[MyNode]    = (MyNode - VertLineId[MyNode]) / NVertLines;
    }
  } else {
    // coordinates mode: we distinguish between vertical line numbering for semi-coarsening and line smoothing
    NumCoords = Ndof / DofsPerNode;

    // reserve temporary memory
    Teuchos::ArrayRCP<LO> TOrigLoc            = Teuchos::arcp<LO>(NumCoords);
    OrigLoc                                   = TOrigLoc.getRawPtr();
    Teuchos::ArrayRCP<coordinate_type> Txtemp = Teuchos::arcp<coordinate_type>(NumCoords);
    xtemp                                     = Txtemp.getRawPtr();
    Teuchos::ArrayRCP<coordinate_type> Tytemp = Teuchos::arcp<coordinate_type>(NumCoords);
    ytemp                                     = Tytemp.getRawPtr();
    Teuchos::ArrayRCP<coordinate_type> Tztemp = Teuchos::arcp<coordinate_type>(NumCoords);
    ztemp                                     = Tztemp.getRawPtr();

    // build vertical line info for semi-coarsening

    // sort coordinates in {x,y,z}vals (returned in {x,y,z}temp) so that we can order things according to lines
    // switch x and y coordinates for semi-coarsening...
    sort_coordinates(NumCoords, OrigLoc, xvals, yvals, zvals, xtemp, ytemp, ztemp, /*true*/ true);

    LO NumBlocks = 0;
    LO index     = 0;

    while (index < NumCoords) {
      xfirst = xtemp[index];
      yfirst = ytemp[index];
      next   = index + 1;
      while ((next != NumCoords) && (xtemp[next] == xfirst) &&
             (ytemp[next] == yfirst))
        next++;
      if (NumBlocks == 0) {
        NumNodesPerVertLine = next - index;
      }
      // The number of vertical lines must be the same on all processors
      // TAW: Sep 14, 2015: or zero as we allow for empty processors.
      // TEUCHOS_TEST_FOR_EXCEPTION(next-index != NumNodesPerVertLine,Exceptions::RuntimeError, "Error code only works for constant block size now!!!\n");
      count = 0;
      for (j = index; j < next; j++) {
        VertLineId[OrigLoc[j]] = NumBlocks;
        LayerId[OrigLoc[j]]    = count++;
      }
      NumBlocks++;
      index = next;
    }
  }

  /* check that everyone was assigned */

  for (i = 0; i < Nnodes; i++) {
    if (VertLineId[i] == -1) {
      GetOStream(Warnings1) << "Warning: did not assign " << i << " to a vertical line?????\n"
                            << std::endl;
    }
    if (LayerId[i] == -1) {
      GetOStream(Warnings1) << "Warning: did not assign " << i << " to a Layer?????\n"
                            << std::endl;
    }
  }

  // TAW: Sep 14 2015: relax plausibility checks as we allow for empty processors
  // MueLu_maxAll(&comm, NumNodesPerVertLine, i);
  // if (NumNodesPerVertLine == -1)  NumNodesPerVertLine = i;
  // TEUCHOS_TEST_FOR_EXCEPTION(NumNodesPerVertLine != i,Exceptions::RuntimeError, "Different processors have different z direction line lengths?\n");

  return NumNodesPerVertLine;
}

/* Private member function to sort coordinates in arrays. This is an expert routine. Do not use or change.*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sort_coordinates(LO numCoords, LO* OrigLoc,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* xvals,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* yvals,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* zvals,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* xtemp,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* ytemp,
                                                                                       typename Teuchos::ScalarTraits<Scalar>::coordinateType* ztemp,
                                                                                       bool flipXY) const {
  if (flipXY == false) {  // for line-smoothing
    for (LO i = 0; i < numCoords; i++) xtemp[i] = xvals[i];
  } else {  // for semi-coarsening
    for (LO i = 0; i < numCoords; i++) xtemp[i] = yvals[i];
  }
  for (LO i = 0; i < numCoords; i++) OrigLoc[i] = i;

  ML_az_dsort2(xtemp, numCoords, OrigLoc);
  if (flipXY == false) {  // for line-smoothing
    for (LO i = 0; i < numCoords; i++) ytemp[i] = yvals[OrigLoc[i]];
  } else {
    for (LO i = 0; i < numCoords; i++) ytemp[i] = xvals[OrigLoc[i]];
  }

  LO index = 0;

  while (index < numCoords) {
    coordinate_type xfirst = xtemp[index];
    LO next                = index + 1;
    while ((next != numCoords) && (xtemp[next] == xfirst))
      next++;
    ML_az_dsort2(&(ytemp[index]), next - index, &(OrigLoc[index]));
    for (LO i = index; i < next; i++) ztemp[i] = zvals[OrigLoc[i]];
    /* One final sort so that the ztemps are in order */
    LO subindex = index;
    while (subindex != next) {
      coordinate_type yfirst = ytemp[subindex];
      LO subnext             = subindex + 1;
      while ((subnext != next) && (ytemp[subnext] == yfirst)) subnext++;
      ML_az_dsort2(&(ztemp[subindex]), subnext - subindex, &(OrigLoc[subindex]));
      subindex = subnext;
    }
    index = next;
  }
}

/* Sort coordinates and additional array accordingly (if provided). This is an expert routine borrowed from ML. Do not change.*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_az_dsort2(typename Teuchos::ScalarTraits<Scalar>::coordinateType dlist[], LocalOrdinal N, LocalOrdinal list2[]) const {
  LO l, r, j, i, flag;
  LO RR2;
  coordinate_type dRR, dK;

  // note: we use that routine for sorting coordinates only. No complex coordinates are assumed...
  typedef Teuchos::ScalarTraits<SC> STS;

  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  dRR = dlist[l - 1];
  dK  = dlist[l - 1];

  if (list2 != NULL) {
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
            if (STS::real(dlist[j]) > STS::real(dlist[j - 1])) j = j + 1;

          if (STS::real(dlist[j - 1]) > STS::real(dK)) {
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
            if (STS::real(dlist[j]) > STS::real(dlist[j - 1])) j = j + 1;
          if (STS::real(dlist[j - 1]) > STS::real(dK)) {
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
}  // namespace MueLu

#endif  // MUELU_LINEDETECTIONFACTORY_DEF_HPP
