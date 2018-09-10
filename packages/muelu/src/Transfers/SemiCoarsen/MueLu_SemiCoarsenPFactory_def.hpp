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
#ifndef MUELU_SEMICOARSENPFACTORY_DEF_HPP
#define MUELU_SEMICOARSENPFACTORY_DEF_HPP

#include <stdlib.h>

#include <Teuchos_LAPACK.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_SemiCoarsenPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("semicoarsen: coarsen rate");
#undef  SET_VALID_ENTRY
    validParamList->set< RCP<const FactoryBase> >("A",               Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",       Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",     Teuchos::null, "Generating factory for coorindates");

    validParamList->set< RCP<const FactoryBase> >("LineDetection_VertLineIds", Teuchos::null, "Generating factory for LineDetection information");
    validParamList->set< RCP<const FactoryBase> >("LineDetection_Layers",      Teuchos::null, "Generating factory for LineDetection information");
    validParamList->set< RCP<const FactoryBase> >("CoarseNumZLayers",          Teuchos::null, "Generating factory for LineDetection information");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");

    Input(fineLevel, "LineDetection_VertLineIds");
    Input(fineLevel, "LineDetection_Layers");
    Input(fineLevel, "CoarseNumZLayers");

    // check whether fine level coordinate information is available.
    // If yes, request the fine level coordinates and generate coarse coordinates
    // during the Build call
    if (fineLevel.GetLevelID() == 0) {
      if (fineLevel.IsAvailable("Coordinates", NoFactory::get())) {
        fineLevel.DeclareInput("Coordinates", NoFactory::get(), this);
        bTransferCoordinates_ = true;
      }
    } else if (bTransferCoordinates_ == true){
      // on coarser levels we check the default factory providing "Coordinates"
      // or the factory declared to provide "Coordinates"
      // first, check which factory is providing coordinate information
      RCP<const FactoryBase> myCoordsFact = GetFactory("Coordinates");
      if (myCoordsFact == Teuchos::null) { myCoordsFact = fineLevel.GetFactoryManager()->GetFactory("Coordinates"); }
      if (fineLevel.IsAvailable("Coordinates", myCoordsFact.get())) {
        fineLevel.DeclareInput("Coordinates", myCoordsFact.get(), this);
        bTransferCoordinates_ = true;
      }
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    // obtain general variables
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");

    // get user-provided coarsening rate parameter (constant over all levels)
    const ParameterList& pL = GetParameterList();
    LO CoarsenRate = as<LO>(pL.get<int>("semicoarsen: coarsen rate"));

    // collect general input data
    LO BlkSize            = A->GetFixedBlockSize();
    RCP<const Map> rowMap = A->getRowMap();
    LO Ndofs              = rowMap->getNodeNumElements();
    LO Nnodes             = Ndofs/BlkSize;

    // collect line detection information generated by the LineDetectionFactory instance
    LO FineNumZLayers = Get< LO >(fineLevel, "CoarseNumZLayers");
    Teuchos::ArrayRCP<LO>     TVertLineId = Get< Teuchos::ArrayRCP<LO> > (fineLevel, "LineDetection_VertLineIds");
    Teuchos::ArrayRCP<LO>     TLayerId    = Get< Teuchos::ArrayRCP<LO> > (fineLevel, "LineDetection_Layers");
    LO* VertLineId = TVertLineId.getRawPtr();
    LO* LayerId    = TLayerId.getRawPtr();

    // generate transfer operator with semicoarsening
    RCP<const Map> theCoarseMap;
    RCP<Matrix>    P;
    RCP<MultiVector> coarseNullspace;
    GO Ncoarse = MakeSemiCoarsenP(Nnodes,FineNumZLayers,CoarsenRate,LayerId,VertLineId,
                               BlkSize, A, P, theCoarseMap, fineNullspace,coarseNullspace);

    // set StridingInformation of P
    if (A->IsView("stridedMaps") == true)
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), theCoarseMap);
    else
      P->CreateView("stridedMaps", P->getRangeMap(), theCoarseMap);

    // Store number of coarse z-layers on the coarse level container
    // This information is used by the LineDetectionAlgorithm
    // TODO get rid of the NoFactory
    
    LO CoarseNumZLayers  =  FineNumZLayers*Ncoarse;
       CoarseNumZLayers  /= Ndofs;
    coarseLevel.Set("NumZLayers", Teuchos::as<LO>(CoarseNumZLayers), MueLu::NoFactory::get());

    // store semicoarsening transfer on coarse level
    Set(coarseLevel, "P", P);

    Set(coarseLevel, "Nullspace", coarseNullspace);

    // transfer coordinates
    if(bTransferCoordinates_) {
      //Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      typedef Xpetra::MultiVector<double,LO,GO,NO> xdMV;
      RCP<xdMV>    fineCoords = Teuchos::null;
      if (fineLevel.GetLevelID() == 0 &&
          fineLevel.IsAvailable("Coordinates", NoFactory::get())) {
        fineCoords = fineLevel.Get< RCP<xdMV> >("Coordinates", NoFactory::get());
      } else {
        RCP<const FactoryBase> myCoordsFact = GetFactory("Coordinates");
        if (myCoordsFact == Teuchos::null) { myCoordsFact = fineLevel.GetFactoryManager()->GetFactory("Coordinates"); }
        if (fineLevel.IsAvailable("Coordinates", myCoordsFact.get())) {
          fineCoords = fineLevel.Get< RCP<xdMV> >("Coordinates", myCoordsFact.get());
        }
      }

      TEUCHOS_TEST_FOR_EXCEPTION(fineCoords==Teuchos::null, Exceptions::RuntimeError, "No Coordinates found provided by the user.");

      TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3, Exceptions::RuntimeError, "Three coordinates arrays must be supplied if line detection orientation not given.");
      ArrayRCP<double> x = fineCoords->getDataNonConst(0);
      ArrayRCP<double> y = fineCoords->getDataNonConst(1);
      ArrayRCP<double> z = fineCoords->getDataNonConst(2);

      // determine the maximum and minimum z coordinate value on the current processor.
      double zval_max = -Teuchos::ScalarTraits<double>::one() / Teuchos::ScalarTraits<double>::sfmin();
      double zval_min =  Teuchos::ScalarTraits<double>::one() / Teuchos::ScalarTraits<double>::sfmin();
      for ( ArrayRCP<double>::iterator it = z.begin(); it != z.end(); ++it) {
        if(*it > zval_max) zval_max = *it;
        if(*it < zval_min) zval_min = *it;
      }

      LO myCoarseZLayers = Teuchos::as<LO>(CoarseNumZLayers);

      ArrayRCP<double> myZLayerCoords = Teuchos::arcp<double>(myCoarseZLayers);
      if(myCoarseZLayers == 1) {
        myZLayerCoords[0] = zval_min;
      } else {
        double dz = (zval_max-zval_min)/(myCoarseZLayers-1);
        for(LO k = 0; k<myCoarseZLayers; ++k) {
          myZLayerCoords[k] = k*dz;
        }
      }

      // Note, that the coarse level node coordinates have to be in vertical ordering according
      // to the numbering of the vertical lines

      // number of vertical lines on current node:
      LO numVertLines = Nnodes / FineNumZLayers;
      LO numLocalCoarseNodes = numVertLines * myCoarseZLayers;

      //std::cout << "rowMap elements: " << rowMap->getNodeNumElements() << std::endl;
      //std::cout << "fineCoords: " << fineCoords->getNodeNumElements() << std::endl;
      //std::cout << "TVertLineId.size(): " << TVertLineId.size() << std::endl;
      //std::cout << "numVertLines=" << numVertLines << std::endl;
      //std::cout << "numLocalCoarseNodes=" << numLocalCoarseNodes << std::endl;

      RCP<const Map>   coarseCoordMap =
          MapFactory::Build (fineCoords->getMap()->lib(),
              Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
              Teuchos::as<size_t>(numLocalCoarseNodes),
              fineCoords->getMap()->getIndexBase(),
              fineCoords->getMap()->getComm());
      RCP<xdMV> coarseCoords   = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(coarseCoordMap, fineCoords->getNumVectors());
      coarseCoords->putScalar(-1.0);
      ArrayRCP<double> cx = coarseCoords->getDataNonConst(0);
      ArrayRCP<double> cy = coarseCoords->getDataNonConst(1);
      ArrayRCP<double> cz = coarseCoords->getDataNonConst(2);

      // loop over all vert line indices (stop as soon as possible)
      LO cntCoarseNodes = 0;
      for( LO vt = 0; vt < TVertLineId.size(); ++vt) {
        //vertical line id in *vt
        LO curVertLineId = TVertLineId[vt];

        if(cx[curVertLineId * myCoarseZLayers] == -1.0 &&
           cy[curVertLineId * myCoarseZLayers] == -1.0) {
          // loop over all local myCoarseZLayers
          for (LO n=0; n<myCoarseZLayers; ++n) {
            cx[curVertLineId * myCoarseZLayers + n] = x[vt];
            cy[curVertLineId * myCoarseZLayers + n] = y[vt];
            cz[curVertLineId * myCoarseZLayers + n] = myZLayerCoords[n];
          }
          cntCoarseNodes += myCoarseZLayers;
        }

        TEUCHOS_TEST_FOR_EXCEPTION(cntCoarseNodes > numLocalCoarseNodes, Exceptions::RuntimeError, "number of coarse nodes is inconsistent.");
        if(cntCoarseNodes == numLocalCoarseNodes) break;
      }

      // set coarse level coordinates
      Set(coarseLevel, "Coordinates", coarseCoords);
    } /* end bool bTransferCoordinates */

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FindCpts(LocalOrdinal const PtsPerLine, LocalOrdinal const CoarsenRate, LocalOrdinal const Thin, LocalOrdinal **LayerCpts) const {
    /*
     * Given the number of points in the z direction (PtsPerLine) and a
     * coarsening rate (CoarsenRate), determine which z-points will serve
     * as Cpts and return the total number of Cpts.
     *
     * Input
     *    PtsPerLine:   Number of fine level points in the z direction
     *
     *    CoarsenRate:  Roughly, number of Cpts  = (PtsPerLine+1)/CoarsenRate - 1
     *
     *    Thin:         Must be either 0 or 1. Thin decides what to do when
     *                  (PtsPerLine+1)/CoarsenRate is not an integer.
     *
     *                    Thin == 0  ==>   ceil() the above fraction
     *                    Thin == 1  ==>   floor() the above fraction
     *
     * Output
     *    LayerCpts     Array where LayerCpts[i] indicates that the
     *                  LayerCpts[i]th fine level layer is a Cpt Layer.
     *                  Note: fine level layers are assumed to be numbered starting
     *                        a one.
     */
    double temp, RestStride, di;
    LO    NCpts, i;
    LO    NCLayers = -1;
    LO    FirstStride;

    temp =  ((double) (PtsPerLine+1))/((double) (CoarsenRate)) - 1.0;
    if  (Thin == 1) NCpts = (LO) ceil(temp);
    else            NCpts = (LO) floor(temp);

    TEUCHOS_TEST_FOR_EXCEPTION(PtsPerLine == 1, Exceptions::RuntimeError, "SemiCoarsenPFactory::FindCpts: cannot coarsen further.");

    if (NCpts < 1) NCpts = 1;

    FirstStride= (LO) ceil( ((double) PtsPerLine+1)/( (double) (NCpts+1)));
    RestStride = ((double) (PtsPerLine-FirstStride+1))/((double) NCpts);

    NCLayers   = (LO) floor((((double) (PtsPerLine-FirstStride+1))/RestStride)+.00001);

    TEUCHOS_TEST_FOR_EXCEPTION(NCLayers != NCpts, Exceptions::RuntimeError, "sizes do not match.");

    di  = (double) FirstStride;
    for (i = 1; i <= NCpts; i++) {
      (*LayerCpts)[i] = (LO) floor(di);
      di += RestStride;
    }

    return(NCLayers);
  }

#define MaxHorNeighborNodes 75

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MakeSemiCoarsenP(LO const Ntotal, LO const nz, LO const CoarsenRate, LO const LayerId[],
                   LO const VertLineId[], LO const DofsPerNode, RCP<Matrix> & Amat, RCP<Matrix>& P,
                   RCP<const Map>& coarseMap, const RCP<MultiVector> fineNullspace,
                   RCP<MultiVector>& coarseNullspace ) const {

    /*
     * Given a CSR matrix (OrigARowPtr, OrigAcols, OrigAvals), information
     * describing the z-layer and vertical line (LayerId and VertLineId)
     * of each matrix block row, a coarsening rate, and dofs/node information,
     * construct a prolongator that coarsening to semicoarsening in the z-direction
     * using something like an operator-dependent grid transfer. In particular,
     * matrix stencils are collapsed to vertical lines. Thus, each vertical line
     * gives rise to a block tridiagonal matrix. BlkRows corresponding to
     * Cpts are replaced by identity matrices. This tridiagonal is solved
     * to determine each interpolation basis functions. Each Blk Rhs corresponds
     * to all zeros except at the corresponding C-pt which has an identity
     *
     * On termination, return the number of local prolongator columns owned by
     * this processor.
     *
     * Note: This code was adapted from a matlab code where offsets/arrays
     *       start from 1. In most parts of the code, this 1 offset is kept
     *       (in some cases wasting the first element of the array). The
     *       input and output matrices of this function has been changed to
     *       have offsets/rows/columns which start from 0. LayerId[] and
     *       VertLineId[] currently start from 1.
     *
     * Input
     * =====
     *    Ntotal       Number of fine level Blk Rows owned by this processor
     *
     *    nz           Number of vertical layers. Note: partitioning must be done
     *                 so that each processor owns an entire vertical line. This
     *                 means that nz is the global number of layers, which should
     *                 be equal to the local number of layers.
     *    CoarsenRate  Rate of z semicoarsening. Smoothed aggregation-like coarsening
     *                 would correspond to CoarsenRate = 3.
     *    LayerId      Array from 0 to Ntotal-1 + Ghost. LayerId(BlkRow) gives the
     *                 layer number associated with the dofs within BlkRow.
     *    VertLineId   Array from 1 to Ntotal, VertLineId(BlkRow) gives a unique
     *                 vertical line id (from 0 to Ntotal/nz-1) of BlkRow. All
     *                 BlkRows associated with nodes along the same vertical line
     *                 in the mesh should have the same LineId.
     *    DofsPerNode  Number of degrees-of-freedom per mesh node.
     *
     *    OrigARowPtr, CSR arrays corresponding to the fine level matrix.
     *    OrigAcols,
     *    OrigAvals
     *
     * Output
     * =====
     *    ParamPptr,   CSR arrays corresponding to the final prolongation matrix.
     *    ParamPcols,
     *    ParamsPvals
     */
    int    NLayers, NVertLines, MaxNnz, NCLayers, MyLine, MyLayer;
    int    *InvLineLayer=NULL, *CptLayers=NULL, StartLayer, NStencilNodes;
    int    BlkRow, dof_j, node_k, *Sub2FullMap=NULL, RowLeng;
    int    i, j, iii, col, count, index, loc, PtRow, PtCol;
    SC     *BandSol=NULL, *BandMat=NULL, TheSum;
    int    *IPIV=NULL, KL, KU, KLU, N, NRHS, LDAB,INFO;
    int    *Pcols;
    size_t *Pptr;
    SC     *Pvals;
    int    MaxStencilSize, MaxNnzPerRow;
    LO     *LayDiff=NULL;
    int    CurRow, LastGuy = -1, NewPtr;
    int    Ndofs;
    int    Nghost;
    LO     *Layerdofs = NULL, *Col2Dof = NULL;

    Teuchos::LAPACK<LO,SC> lapack;

    char notrans[3];
    notrans[0] = 'N';
    notrans[1] = 'N';


    MaxNnzPerRow = MaxHorNeighborNodes*DofsPerNode*3;
    Teuchos::ArrayRCP<LO> TLayDiff = Teuchos::arcp<LO>(1+MaxNnzPerRow); LayDiff = TLayDiff.getRawPtr();

    Nghost = Amat->getColMap()->getNodeNumElements() - Amat->getDomainMap()->getNodeNumElements();
    if (Nghost < 0) Nghost = 0;
    Teuchos::ArrayRCP<LO> TLayerdofs= Teuchos::arcp<LO>(Ntotal*DofsPerNode+Nghost+1); Layerdofs = TLayerdofs.getRawPtr();
    Teuchos::ArrayRCP<LO> TCol2Dof= Teuchos::arcp<LO>(Ntotal*DofsPerNode+Nghost+1); Col2Dof= TCol2Dof.getRawPtr();

    RCP<Xpetra::Vector<LO,LO,GO,NO> > localdtemp = Xpetra::VectorFactory<LO,LO,GO,NO>::Build(Amat->getDomainMap());
    RCP<Xpetra::Vector<LO,LO,GO,NO> > dtemp      = Xpetra::VectorFactory<LO,LO,GO,NO>::Build(Amat->getColMap());
    ArrayRCP<LO> valptr= localdtemp->getDataNonConst(0);

    for (i = 0; i < Ntotal*DofsPerNode; i++)
      valptr[i]= LayerId[i/DofsPerNode];

    RCP< const Import> importer;
    importer = Amat->getCrsGraph()->getImporter();
    if (importer == Teuchos::null) {
      importer = ImportFactory::Build(Amat->getDomainMap(), Amat->getColMap());
    }
    dtemp->doImport(*localdtemp, *(importer), Xpetra::INSERT);

    valptr= dtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Layerdofs[i]= valptr[i];
    valptr= localdtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode;        i++) valptr[i]= i%DofsPerNode;
    dtemp->doImport(*localdtemp, *(importer), Xpetra::INSERT);
    valptr= dtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Col2Dof[i]= valptr[i];

    if (Ntotal != 0) {
      NLayers   = LayerId[0];
      NVertLines= VertLineId[0];
    }
    else { NLayers = -1; NVertLines = -1; }

    for (i = 1; i < Ntotal; i++) {
      if ( VertLineId[i] > NVertLines ) NVertLines = VertLineId[i];
      if ( LayerId[i]    >   NLayers  ) NLayers    = LayerId[i];
    }
    NLayers++;
    NVertLines++;

    /*
     * Make an inverse map so that we can quickly find the dof
     * associated with a particular vertical line and layer.
     */

    Teuchos::ArrayRCP<LO> TInvLineLayer= Teuchos::arcp<LO>(1+NVertLines*NLayers); InvLineLayer = TInvLineLayer.getRawPtr();
    for (i=0; i < Ntotal; i++) {
      InvLineLayer[ VertLineId[i]+1+LayerId[i]*NVertLines ] = i;
    }

    /*
     * Determine coarse layers where injection will be applied.
     */

    Teuchos::ArrayRCP<LO> TCptLayers = Teuchos::arcp<LO>(nz+1); CptLayers = TCptLayers.getRawPtr();
    NCLayers = FindCpts(nz,CoarsenRate,0, &CptLayers);

    /*
     * Compute the largest possible interpolation stencil width based
     * on the location of the Clayers. This stencil width is actually
     * nodal (i.e. assuming 1 dof/node). To get the true max stencil width
     * one needs to multiply this by DofsPerNode.
     */

    if  (NCLayers < 2) MaxStencilSize = nz;
    else MaxStencilSize = CptLayers[2];

    for (i = 3; i <= NCLayers; i++) {
      if (MaxStencilSize < CptLayers[i]- CptLayers[i-2])
        MaxStencilSize = CptLayers[i]- CptLayers[i-2];
    }
    if (NCLayers > 1) {
      if (MaxStencilSize < nz - CptLayers[NCLayers-1]+1)
        MaxStencilSize =  nz - CptLayers[NCLayers-1]+1;
    }

    /*
     * Allocate storage associated with solving a banded sub-matrix needed to
     * determine the interpolation stencil. Note: we compute interpolation stencils
     * for all dofs within a node at the same time, and so the banded solution
     * must be large enough to hold all DofsPerNode simultaneously.
     */

    Teuchos::ArrayRCP<LO> TSub2FullMap= Teuchos::arcp<LO>((MaxStencilSize+1)*DofsPerNode); Sub2FullMap= TSub2FullMap.getRawPtr();
    Teuchos::ArrayRCP<SC> TBandSol= Teuchos::arcp<SC>((MaxStencilSize+1)*DofsPerNode*DofsPerNode); BandSol = TBandSol.getRawPtr();
    /*
     * Lapack variables. See comments for dgbsv().
     */
    KL     = 2*DofsPerNode-1;
    KU     = 2*DofsPerNode-1;
    KLU    = KL+KU;
    LDAB   = 2*KL+KU+1;
    NRHS   = DofsPerNode;
    Teuchos::ArrayRCP<SC> TBandMat= Teuchos::arcp<SC>(LDAB*MaxStencilSize*DofsPerNode+1); BandMat = TBandMat.getRawPtr();
    Teuchos::ArrayRCP<LO> TIPIV= Teuchos::arcp<LO>((MaxStencilSize+1)*DofsPerNode); IPIV = TIPIV.getRawPtr();

    /*
     * Allocate storage for the final interpolation matrix. Note: each prolongator
     * row might have entries corresponding to at most two nodes.
     * Note: the total fine level dofs equals DofsPerNode*Ntotal and the max
     *       nnz per prolongator row is DofsPerNode*2.
     */

    Ndofs  = DofsPerNode*Ntotal;
    MaxNnz = 2*DofsPerNode*Ndofs;

    RCP<const Map> rowMap = Amat->getRowMap();
    int GNdofs= rowMap->getGlobalNumElements();

    std::vector<size_t> stridingInfo_;
    stridingInfo_.push_back(DofsPerNode);

    coarseMap = StridedMapFactory::Build(rowMap->lib(),
        (NCLayers*GNdofs)/nz,
        NCLayers*NVertLines*DofsPerNode,
        0, /* index base */
        stridingInfo_,
        rowMap->getComm(),
        -1, /* strided block id */
        0, /* domain gid offset */
        rowMap->getNode());


    //coarseMap = MapFactory::createContigMapWithNode(rowMap->lib(),(NCLayers*GNdofs)/nz, NCLayers*NVertLines*DofsPerNode,(rowMap->getComm()), rowMap->getNode());


    P       = rcp(new CrsMatrixWrap(rowMap, coarseMap , 0, Xpetra::StaticProfile));
    coarseNullspace = MultiVectorFactory::Build(coarseMap, fineNullspace->getNumVectors());


    Teuchos::ArrayRCP<SC> TPvals= Teuchos::arcp<SC>(1+MaxNnz);               Pvals= TPvals.getRawPtr();
    Teuchos::ArrayRCP<size_t> TPptr = Teuchos::arcp<size_t>(DofsPerNode*(2+Ntotal)); Pptr = TPptr.getRawPtr();
    Teuchos::ArrayRCP<LO>     TPcols= Teuchos::arcp<LO>(1+MaxNnz);                   Pcols= TPcols.getRawPtr();

    Pptr[0] = 0; Pptr[1] = 0;

    TEUCHOS_TEST_FOR_EXCEPTION(Pcols == NULL, Exceptions::RuntimeError, "MakeSemiCoarsenP: Not enough space \n");

    /*
     * Setup P's rowptr as if each row had its maximum of 2*DofsPerNode nonzeros.
     * This will be useful while filling up P, and then later we will squeeze out
     * the unused nonzeros locations.
     */

    for (i = 1; i <= MaxNnz; i++) Pcols[i] = -1;  /* mark all entries as unused */
    count = 1;
    for (i = 1; i <= DofsPerNode*Ntotal+1; i++) {
      Pptr[i]  = count;
      count   += (2*DofsPerNode);
    }

    /*
     * Build P column by column. The 1st block column corresponds to the 1st coarse
     * layer and the first line. The 2nd block column corresponds to the 2nd coarse
     * layer and the first line. The NCLayers+1 block column corresponds to the
     * 1st coarse layer and the 2nd line, etc.
     */


    col = 0;
    for (MyLine=1; MyLine <= NVertLines; MyLine += 1) {
      for (iii=1; iii <= NCLayers;  iii+= 1) {
        col = col+1;
        MyLayer = CptLayers[iii];

        /*
         * StartLayer gives the layer number of the lowest layer that
         * is nonzero in the interpolation stencil that is currently
         * being computed. Normally, if we are not near a boundary this
         * is simply CptsLayers[iii-1]+1.
         *
         * NStencilNodes indicates the number of nonzero nodes in the
         * interpolation stencil that is currently being computed. Normally,
         * if we are not near a boundary this is CptLayers[iii+1]-StartLayer.
         */

        if (iii !=    1    ) StartLayer = CptLayers[iii-1]+1;
        else                 StartLayer = 1;

        if (iii != NCLayers) NStencilNodes= CptLayers[iii+1]-StartLayer;
        else                 NStencilNodes= NLayers - StartLayer + 1;


        N = NStencilNodes*DofsPerNode;

        /*
         *  dgbsv() does not require that the first KL rows be initialized,
         *  so we could avoid zeroing out some entries?
         */

        for (i = 0; i < NStencilNodes*DofsPerNode*DofsPerNode; i++)
          BandSol[ i] = 0.0;
        for (i = 0; i < LDAB*N; i++) BandMat[ i] = 0.0;

        /*
         *  Fill BandMat and BandSol (which is initially the rhs) for each
         *  node in the interpolation stencil that is being computed.
         */

        for (node_k=1; node_k <= NStencilNodes ; node_k++) {

          /*  Map a Line and Layer number to a BlkRow in the fine level  matrix
           *  and record the mapping from the sub-system to the BlkRow of the
           *  fine level matrix.
           */
          BlkRow  = InvLineLayer[MyLine+(StartLayer+node_k-2)*NVertLines]+1;
          Sub2FullMap[node_k] = BlkRow;

          /* Two cases:
           *    1) the current layer is not a Cpoint layer. In this case we
           *       want to basically stick the matrix couplings to other
           *       nonzero stencil rows into the band matrix. One way to do
           *       this is to include couplings associated with only MyLine
           *       and ignore all the other couplings. However, what we do
           *       instead is to sum all the coupling at each layer participating
           *       in this interpolation stencil and stick this sum into BandMat.
           *    2) the current layer is a Cpoint layer and so we
           *       stick an identity block in BandMat and rhs.
           */
          if (StartLayer+node_k-1 != MyLayer) {
            for (int dof_i=0; dof_i < DofsPerNode; dof_i++) {

              j = (BlkRow-1)*DofsPerNode+dof_i;
              ArrayView<const LO> AAcols;
              ArrayView<const SC> AAvals;
              Amat->getLocalRowView(j, AAcols, AAvals);
              const int *Acols    = AAcols.getRawPtr();
              const SC *Avals = AAvals.getRawPtr();
              RowLeng = AAvals.size();

              TEUCHOS_TEST_FOR_EXCEPTION(RowLeng >= MaxNnzPerRow, Exceptions::RuntimeError, "MakeSemiCoarsenP: recompile with larger Max(HorNeighborNodes)\n");

              for (i = 0; i < RowLeng; i++) {
                LayDiff[i]  = Layerdofs[Acols[i]]-StartLayer-node_k+2;

                /* This is the main spot where there might be off- */
                /* processor communication. That is, when we       */
                /* average the stencil in the horizontal direction,*/
                /* we need to know the layer id of some            */
                /* neighbors that might reside off-processor.      */
              }
              PtRow = (node_k-1)*DofsPerNode+dof_i+1;
              for (dof_j=0; dof_j < DofsPerNode; dof_j++) {
                PtCol = (node_k-1)*DofsPerNode+dof_j + 1;
                /* Stick in entry corresponding to Mat(PtRow,PtCol) */
                /* see dgbsv() comments for matrix format.          */
                TheSum = 0.0;
                for (i = 0; i < RowLeng; i++) {
                  if ((LayDiff[i] == 0)  && (Col2Dof[Acols[i]] == dof_j))
                    TheSum += Avals[i];
                }
                index = LDAB*(PtCol-1)+KLU+PtRow-PtCol;
                BandMat[index] = TheSum;
                if (node_k != NStencilNodes) {
                  /* Stick Mat(PtRow,PtCol+DofsPerNode) entry  */
                  /* see dgbsv() comments for matrix format.  */
                  TheSum = 0.0;
                  for (i = 0; i < RowLeng; i++) {
                    if ((LayDiff[i] == 1) &&(Col2Dof[Acols[i]]== dof_j))
                      TheSum += Avals[i];
                  }
                  j = PtCol+DofsPerNode;
                  index=LDAB*(j-1)+KLU+PtRow-j;
                  BandMat[index] = TheSum;
                }
                if (node_k != 1) {
                  /* Stick Mat(PtRow,PtCol-DofsPerNode) entry  */
                  /* see dgbsv() comments for matrix format.  */
                  TheSum = 0.0;
                  for (i = 0; i < RowLeng; i++) {
                    if ((LayDiff[i]== -1) &&(Col2Dof[Acols[i]]== dof_j))
                      TheSum += Avals[i];
                  }
                  j = PtCol-DofsPerNode;
                  index=LDAB*(j-1)+KLU+PtRow-j;
                  BandMat[index] = TheSum;
                }
              }
            }
          }
          else {

             /* inject the null space */
            // int FineStride  = Ntotal*DofsPerNode;
            // int CoarseStride= NVertLines*NCLayers*DofsPerNode;
            for (int k = 0; k < fineNullspace->getNumVectors(); k++) {
              Teuchos::ArrayRCP<SC> OneCNull = coarseNullspace->getDataNonConst(k);
              Teuchos::ArrayRCP<SC> OneFNull = fineNullspace->getDataNonConst(k);
              for (int dof_i = 0; dof_i < DofsPerNode; dof_i++) {
                 OneCNull[(   col-1)*DofsPerNode+dof_i] = OneFNull[ (BlkRow-1)*DofsPerNode+dof_i];
              }
            }

            for (int dof_i = 0; dof_i < DofsPerNode; dof_i++) {
              /* Stick Mat(PtRow,PtRow) and Rhs(PtRow,dof_i+1) */
              /* see dgbsv() comments for matrix format.     */
              PtRow = (node_k-1)*DofsPerNode+dof_i+1;
              index = LDAB*(PtRow-1)+KLU;
              BandMat[index] = 1.0;
              BandSol[(dof_i)*DofsPerNode*NStencilNodes+PtRow-1] = 1.;
            }
          }
        }

        /* Solve banded system and then stick result in Pmatrix arrays */

        lapack.GBTRF( N, N, KL, KU, BandMat, LDAB, IPIV, &INFO);

        TEUCHOS_TEST_FOR_EXCEPTION(INFO != 0, Exceptions::RuntimeError, "Lapack band factorization failed");

        lapack.GBTRS(notrans[0], N, KL, KU, NRHS, BandMat, LDAB, IPIV,
                     BandSol, N, &INFO );

        TEUCHOS_TEST_FOR_EXCEPTION(INFO != 0, Exceptions::RuntimeError, "Lapack band solve back substitution failed");

        for (dof_j=0; dof_j < DofsPerNode; dof_j++) {
          for (int dof_i=0; dof_i < DofsPerNode; dof_i++) {
            for (i =1; i <= NStencilNodes ; i++) {
              index = (Sub2FullMap[i]-1)*DofsPerNode+dof_i+1;
              loc = Pptr[index];
              Pcols[loc] = (col-1)*DofsPerNode+dof_j+1;
              Pvals[loc] = BandSol[dof_j*DofsPerNode*NStencilNodes +
                  (i-1)*DofsPerNode + dof_i];
              Pptr[index]= Pptr[index] + 1;
            }
          }
        }
      }
    }

    /*
     * Squeeze the -1's out of the columns. At the same time convert Pcols
     * so that now the first column is numbered '0' as opposed to '1'.
     * Also, the arrays Pcols and Pvals should now use the zeroth element
     * as opposed to just starting with the first element. Pptr will be
     * fixed in the for loop below so that Pptr[0] = 0, etc.
     */
    CurRow = 1;
    NewPtr = 1;
    for (size_t ii=1; ii <= Pptr[Ntotal*DofsPerNode]-1; ii++) {
      if (ii == Pptr[CurRow]) {
        Pptr[CurRow] = LastGuy;
        CurRow++;
        while (ii > Pptr[CurRow]) {
          Pptr[CurRow] = LastGuy;
          CurRow++;
        }
      }
      if (Pcols[ii] != -1) {
        Pcols[NewPtr-1] = Pcols[ii]-1;   /* these -1's fix the offset and */
        Pvals[NewPtr-1] = Pvals[ii];     /* start using the zeroth element */
        LastGuy = NewPtr;
        NewPtr++;
      }
    }
    for (i = CurRow; i <= Ntotal*DofsPerNode; i++) Pptr[CurRow] = LastGuy;

    /* Now move the pointers so that they now point to the beginning of each
     * row as opposed to the end of each row
     */
    for (i=-Ntotal*DofsPerNode+1; i>= 2 ; i--) {
      Pptr[i-1] = Pptr[i-2];  /* extra -1 added to start from 0 */
    }
    Pptr[0] = 0;

    ArrayRCP<size_t>  rcpRowPtr;
    ArrayRCP<LO>      rcpColumns;
    ArrayRCP<SC>      rcpValues;

    RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    LO nnz =  Pptr[Ndofs];
    PCrs->allocateAllValues(nnz, rcpRowPtr, rcpColumns, rcpValues);

    ArrayView<size_t> rowPtr  = rcpRowPtr();
    ArrayView<LO>     columns = rcpColumns();
    ArrayView<SC>     values  = rcpValues();

    // copy data over

    for (LO ii = 0; ii <= Ndofs; ii++) rowPtr[ii] = Pptr[ii];
    for (LO ii = 0; ii < nnz; ii++)    columns[ii] = Pcols[ii];
    for (LO ii = 0; ii < nnz; ii++) values[ii]  = Pvals[ii];
    PCrs->setAllValues(rcpRowPtr, rcpColumns, rcpValues);
    PCrs->expertStaticFillComplete(coarseMap, Amat->getDomainMap());


    return NCLayers*NVertLines*DofsPerNode;
  }
} //namespace MueLu

#define MUELU_SEMICOARSENPFACTORY_SHORT
#endif // MUELU_SEMICOARSENPFACTORY_DEF_HPP
