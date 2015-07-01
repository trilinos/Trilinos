// @HEADER
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
#ifndef MUELU_LINEDETECTIONFACTORY_DEF_HPP
#define MUELU_LINEDETECTIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_LineDetectionFactory_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",               Teuchos::null, "Generating factory of the matrix A");
    //validParamList->set< RCP<const FactoryBase> >("Nullspace",       Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",     Teuchos::null, "Generating factory for coorindates");

    validParamList->set< std::string > ("linedetection: orientation", "vertical", "Line orientation: can be either 'vertical', 'horizontal' or 'coordinates'");
    validParamList->set< LO > ("linedetection: num layers", 10, "Line detection: number of layers on finest level. Alternatively, set the number of layers on the finest level as \"NumZLayers\" in the finest level container class.");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    // what about Coordinates

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

    LO               NumZDir = 0, Zorientation = GRID_SUPPLIED;
    RCP<MultiVector> fineCoords;
    ArrayRCP<Scalar> x, y, z;
    Scalar           *xptr = NULL, *yptr = NULL, *zptr = NULL;

    // obtain general variables
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (currentLevel, "A");

    const ParameterList& pL = GetParameterList();
    const std::string lineOrientation = pL.get<std::string>("linedetection: orientation");
    GetOStream(Runtime1) << "Line detection mode = " << lineOrientation << std::endl;

    // interpret "line orientation" parameter provided by the user on the finest level
    if(currentLevel.GetLevelID() == 0) {
      if(lineOrientation=="vertical")
        Zorientation = VERTICAL;
      else if (lineOrientation=="horizontal")
        Zorientation = HORIZONTAL;
      else if (lineOrientation=="coordinates")
        Zorientation = GRID_SUPPLIED;
      else
        TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "LineDetectionFactory: The parameter 'semicoarsen: line orientation' must be either 'vertical', 'horizontal' or 'coordinates'.");
    } else {
      // on coarse levels the line orientation is internally fixed to be vertical
      Zorientation = VERTICAL;
    }

    // obtain number of z layers (variable over levels)
    // This information is user-provided on the finest level and transferred to the coarser
    // levels by the SemiCoarsenPFactor using the internal "NumZLayers" variable.
    if(currentLevel.GetLevelID() == 0) {
      if(currentLevel.IsAvailable("NumZLayers", NoFactory::get())) {
        NumZDir = currentLevel.Get<LO>("NumZLayers", NoFactory::get()); //obtain info
        GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << " (information from Level(0))" << std::endl;
      } else {
        NumZDir = pL.get<LO>("linedetection: num layers");
        GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << " (information provided by user through 'line detection: num layers')" << std::endl;
      }
    } else {
      // TODO get rid of NoFactory here and use SemiCoarsenPFactory as source of NumZLayers instead.
      if(currentLevel.IsAvailable("NumZLayers", NoFactory::get())) {
        NumZDir = currentLevel.Get<LO>("NumZLayers", NoFactory::get()); //obtain info
        GetOStream(Runtime1) << "Number of layers for line detection: " << NumZDir << std::endl;
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "LineDetectionFactory: BuildP: No NumZLayers variable found. This cannot be.");
      }
    }

    GetOStream(Runtime1) << "Found " << NumZDir << " layers" << std::endl;

    // plausibility check and further variable collection
    if (Zorientation == GRID_SUPPLIED) { // On finest level, fetch user-provided coordinates if available...
      // NOTE: this code ist not tested!
      bool CoordsAvail = currentLevel.IsAvailable("Coordinates");

      if (CoordsAvail == false) {
        if (currentLevel.GetLevelID() == 0)
          throw Exceptions::RuntimeError("Coordinates must be supplied if line detection orientation not given.");
        else
          throw Exceptions::RuntimeError("Coordinates not generated by previous invocation of LineDetectionFactory's BuildP() method.");
      }
      fineCoords = Get< RCP<MultiVector> > (currentLevel, "Coordinates");
      TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3, Exceptions::RuntimeError, "Three coordinates arrays must be supplied if line detection orientation not given.");
      x = fineCoords->getDataNonConst(0);
      y = fineCoords->getDataNonConst(1);
      z = fineCoords->getDataNonConst(2);
      xptr = x.getRawPtr();
      yptr = y.getRawPtr();
      zptr = z.getRawPtr();
    }

    // collect common information
    LO BlkSize = A->GetFixedBlockSize();
    //TEUCHOS_TEST_FOR_EXCEPTION(BlkSize != 1, Exceptions::RuntimeError, "Block size > 1 has not been implemented");

    RCP<const Map> rowMap = A->getRowMap();
    LO Ndofs   = rowMap->getNodeNumElements();
    LO Nnodes  = Ndofs/BlkSize;

    // perform line detection
    if (NumZDir > 0) {
      LO   *LayerId, *VertLineId;
      //Teuchos::ArrayRCP<LO>     TLayerId   = Teuchos::arcp<LO>(Nnodes+1);  LayerId   = TLayerId.getRawPtr();
      //Teuchos::ArrayRCP<LO>     TVertLineId= Teuchos::arcp<LO>(Nnodes+1);  VertLineId= TVertLineId.getRawPtr();
      Teuchos::ArrayRCP<LO>     TLayerId   = Teuchos::arcp<LO>(Nnodes);  LayerId   = TLayerId.getRawPtr();
      Teuchos::ArrayRCP<LO>     TVertLineId= Teuchos::arcp<LO>(Nnodes);  VertLineId= TVertLineId.getRawPtr();

      NumZDir = ML_compute_line_info(LayerId, VertLineId, Ndofs, BlkSize,
                                     Zorientation, NumZDir,xptr,yptr,zptr, *(rowMap->getComm()));
      //it is NumZDir=NCLayers*NVertLines*DofsPerNode;

      // store output data on current level
      // The line detection data is used by the SemiCoarsenPFactory and the line smoothers in Ifpack/Ifpack2
      Set(currentLevel, "CoarseNumZLayers", NumZDir);
      Set(currentLevel, "LineDetection_Layers", TLayerId);
      Set(currentLevel, "LineDetection_VertLineIds", TVertLineId);
    } else {
      Teuchos::ArrayRCP<LO>     TLayerId   = Teuchos::arcp<LO>(0);
      Teuchos::ArrayRCP<LO>     TVertLineId= Teuchos::arcp<LO>(0);

      // store output data on current level
      // The line detection data is used by the SemiCoarsenPFactory and the line smoothers in Ifpack/Ifpack2
      Set(currentLevel, "CoarseNumZLayers", NumZDir);
      Set(currentLevel, "LineDetection_Layers", TLayerId);
      Set(currentLevel, "LineDetection_VertLineIds", TVertLineId);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_compute_line_info(LocalOrdinal LayerId[], LocalOrdinal VertLineId[], LocalOrdinal Ndof, LocalOrdinal DofsPerNode, LocalOrdinal MeshNumbering, LocalOrdinal NumNodesPerVertLine, Scalar *xvals, Scalar *yvals, Scalar *zvals, const Teuchos::Comm<int>& comm) const {

    LO    Nnodes, NVertLines, MyNode;
    LO    NumCoords, NumBlocks, index, next, subindex, subnext;
    SC xfirst, yfirst;
    SC *xtemp, *ytemp, *ztemp;
    LO    *OrigLoc;
    LO    i,j,count;
    LO    RetVal;
    //LO    mypid; // Not used

    //mypid = comm.getRank();
    RetVal = 0;
    if ((MeshNumbering != VERTICAL) && (MeshNumbering != HORIZONTAL)) {
      if ( (xvals == NULL) || (yvals == NULL) || (zvals == NULL)) RetVal = -1;
    }
    else {
      if  (NumNodesPerVertLine == -1)                     RetVal = -4;
      if ( ((Ndof/DofsPerNode)%NumNodesPerVertLine) != 0) RetVal = -3;
    }
    if ( (Ndof%DofsPerNode) != 0) RetVal = -2;

    TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -1, Exceptions::RuntimeError, "Not semicoarsening as no mesh numbering information or coordinates are given\n");
    TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -4, Exceptions::RuntimeError, "Not semicoarsening as the number of z nodes is not given.\n");
    TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -3, Exceptions::RuntimeError, "Not semicoarsening as the total number of nodes is not evenly divisible by the number of z direction nodes .\n");
    TEUCHOS_TEST_FOR_EXCEPTION(RetVal == -2, Exceptions::RuntimeError, "Not semicoarsening as something is off with the number of degrees-of-freedom per node.\n");

    Nnodes = Ndof/DofsPerNode;

    for (MyNode = 0; MyNode < Nnodes;  MyNode++) VertLineId[MyNode]= -1;
    for (MyNode = 0; MyNode < Nnodes;  MyNode++) LayerId[MyNode]   = -1;


    if (MeshNumbering == VERTICAL) {
      for (MyNode = 0; MyNode < Nnodes; MyNode++) {
        LayerId[MyNode]= MyNode%NumNodesPerVertLine;
        VertLineId[MyNode]= (MyNode- LayerId[MyNode])/NumNodesPerVertLine;
      }
    }
    else if (MeshNumbering == HORIZONTAL) {
      NVertLines = Nnodes/NumNodesPerVertLine;
      for (MyNode = 0; MyNode < Nnodes; MyNode++) {
        VertLineId[MyNode]   = MyNode%NVertLines;
        LayerId[MyNode]   = (MyNode- VertLineId[MyNode])/NVertLines;
      }
    }
    else {


      NumCoords = Ndof/DofsPerNode;

      /* sort coordinates so that we can order things according to lines */

      Teuchos::ArrayRCP<LO> TOrigLoc= Teuchos::arcp<LO>(NumCoords+1);       OrigLoc= TOrigLoc.getRawPtr();
      Teuchos::ArrayRCP<SC> Txtemp  = Teuchos::arcp<SC>(NumCoords+1);       xtemp  = Txtemp.getRawPtr();
      Teuchos::ArrayRCP<SC> Tytemp  = Teuchos::arcp<SC>(NumCoords+1);       ytemp  = Tytemp.getRawPtr();
      Teuchos::ArrayRCP<SC> Tztemp  = Teuchos::arcp<SC>(NumCoords+1);       ztemp  = Tztemp.getRawPtr();

      TEUCHOS_TEST_FOR_EXCEPTION(ztemp == NULL, Exceptions::RuntimeError, "Not enough memory for line algorithms");
      for (i = 0; i < NumCoords; i++) ytemp[i]= yvals[i];
      for (i = 0; i < NumCoords; i++) OrigLoc[i]= i;

      ML_az_dsort2(ytemp,NumCoords,OrigLoc);
      for (i = 0; i < NumCoords; i++) xtemp[i]= xvals[OrigLoc[i]];

      index = 0;

      while ( index < NumCoords ) {
        yfirst = ytemp[index];
        next   = index+1;
        while ( (next != NumCoords) && (ytemp[next] == yfirst))
          next++;
        ML_az_dsort2(&(xtemp[index]),next-index,&(OrigLoc[index]));
        for (i = index; i < next; i++) ztemp[i]= zvals[OrigLoc[i]];
        /* One final sort so that the ztemps are in order */
        subindex = index;
        while (subindex != next) {
          xfirst = xtemp[subindex]; subnext = subindex+1;
          while ( (subnext != next) && (xtemp[subnext] == xfirst)) subnext++;
          ML_az_dsort2(&(ztemp[subindex]),subnext-subindex,&(OrigLoc[subindex]));
          subindex = subnext;
        }
        index = next;
      }

      /* go through each vertical line and populate blockIndices so all   */
      /* dofs within a PDE within a vertical line correspond to one block.*/

      NumBlocks = 0;
      index = 0;

      while ( index < NumCoords ) {
        xfirst = xtemp[index];  yfirst = ytemp[index];
        next = index+1;
        while ( (next != NumCoords) && (xtemp[next] == xfirst) &&
                (ytemp[next] == yfirst))
          next++;
        if (NumBlocks == 0) NumNodesPerVertLine = next-index;
        TEUCHOS_TEST_FOR_EXCEPTION(next-index != NumNodesPerVertLine,Exceptions::RuntimeError, "Error code only works for constant block size now!!!\n");
        count = 0;
        for (j= index; j < next; j++) {
          VertLineId[OrigLoc[j]] = NumBlocks;
          LayerId[OrigLoc[j]] = count++;
        }
        NumBlocks++;
        index = next;
      }
    }

    /* check that everyone was assigned */

    for (i = 0; i < Nnodes;  i++) {
      if (VertLineId[i] == -1) {
        std::cout << "Warning: did not assign " << i << " to a vertical line?????\n" << std::endl;
      }
      if (LayerId[i] == -1) {
        std::cout << "Warning: did not assign " << i << " to a Layer?????\n" << std::endl;
      }
    }
    MueLu_maxAll(&comm, NumNodesPerVertLine, i);
    if (NumNodesPerVertLine == -1)  NumNodesPerVertLine = i;

    TEUCHOS_TEST_FOR_EXCEPTION(NumNodesPerVertLine != i,Exceptions::RuntimeError, "Different processors have different z direction line lengths?\n");
    return NumNodesPerVertLine;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_az_dsort2(Scalar dlist[], LocalOrdinal N, LocalOrdinal list2[]) const {
    LO l, r, j, i, flag;
    LO RR2;
    SC       dRR, dK;
    typedef Teuchos::ScalarTraits<SC> STS;

    if (N <= 1) return;

    l    = N / 2 + 1;
    r    = N - 1;
    l    = l - 1;
    dRR  = dlist[l - 1];
    dK   = dlist[l - 1];

    if (list2 != NULL) {
      RR2 = list2[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (STS::magnitude(dlist[j]) > STS::magnitude(dlist[j - 1])) j = j + 1;

            if (STS::magnitude(dlist[j - 1]) > STS::magnitude(dK)) {
              dlist[ i - 1] = dlist[ j - 1];
              list2[i - 1] = list2[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }
        dlist[ i - 1] = dRR;
        list2[i - 1] = RR2;

        if (l == 1) {
          dRR  = dlist [r];
          RR2 = list2[r];
          dK = dlist[r];
          dlist[r ] = dlist[0];
          list2[r] = list2[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          dRR  = dlist[ l - 1];
          RR2 = list2[l - 1];
          dK   = dlist[l - 1];
        }
      }
      dlist[ 0] = dRR;
      list2[0] = RR2;
    }
    else {
      while (r != 0) {
        j = l;
        flag = 1;
        while (flag == 1) {
          i = j;
          j = j + j;
          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (STS::magnitude(dlist[j]) > STS::magnitude(dlist[j - 1])) j = j + 1;
            if (STS::magnitude(dlist[j - 1]) > STS::magnitude(dK)) {
              dlist[ i - 1] = dlist[ j - 1];
            }
            else {
              flag = 0;
            }
          }
        }
        dlist[ i - 1] = dRR;
        if (l == 1) {
          dRR  = dlist [r];
          dK = dlist[r];
          dlist[r ] = dlist[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          dRR  = dlist[ l - 1];
          dK   = dlist[l - 1];
        }
      }
      dlist[ 0] = dRR;
    }

  }
} //namespace MueLu

#endif // MUELU_LINEDETECTIONFACTORY_DEF_HPP
