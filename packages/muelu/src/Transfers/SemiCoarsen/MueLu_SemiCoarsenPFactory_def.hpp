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

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_SemiCoarsenPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"
#include <Teuchos_LAPACK.hpp>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("semicoarsen: coarsen rate");
#undef  SET_VALID_ENTRY
    validParamList->set<ArrayRCP<LO> >("SemiCoarsenInfo",  Teuchos::null, "Generating factory of the array SemiCoarsenInfo");

    validParamList->set< RCP<const FactoryBase> >("A",            Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",    Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",  Teuchos::null, "Generating factory for coorindates");


    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    //    Input(fineLevel, "Coordinates");   // rst: this didn't work last time I tried?
    //    Input(fineLevel, "SemiCoarsenInfo");    // rst: this didn't work last time I tried?

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    LO               NumZDir, Zorientation;
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    ArrayRCP<LO>     FineSemiInfo = fineLevel.Get<ArrayRCP<LO> >("SemiCoarsenInfo"); // rst: why do I need this style of Get?
    //    ArrayRCP<LO>     FineSemiInfo = Get< ArrayRCP<LO> >(fineLevel,"SemiCoarsenInfo"); // rst: why do I need this style of Get?

    NumZDir      = FineSemiInfo[NUM_ZPTS];
    Zorientation = FineSemiInfo[ORIENTATION];


    LO BlkSize = A->GetFixedBlockSize();
    TEUCHOS_TEST_FOR_EXCEPTION(BlkSize != 1, Exceptions::RuntimeError, "Block size > 1 has not been implemented");

    RCP<const Map> rowMap = A->getRowMap();

    LO   Ndofs, Nnodes, *LayerId, *VertLineId, CoarsenRate;
    GO   Ncoarse;

    const ParameterList &pL = GetParameterList();
    CoarsenRate = as<LO>(pL.get<int>("semicoarsen: coarsen rate"));

    Ndofs   = rowMap->getNodeNumElements();
    Nnodes  = Ndofs/BlkSize;


    RCP<MultiVector> fineCoords;
    ArrayRCP<Scalar> x, y, z;
    SC *xptr = NULL, *yptr = NULL, *zptr = NULL;

    if ( (Zorientation != VERTICAL) &&  (Zorientation != HORIZONTAL) ) {
      TEUCHOS_TEST_FOR_EXCEPTION(!fineLevel.IsAvailable("Coordinates"), Exceptions::RuntimeError, "Coordinates must be supplied if semicoarsen orientation not given.");
      fineCoords = Get< RCP<MultiVector> > (fineLevel, "Coordinates");
      TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3, Exceptions::RuntimeError, "Three coordinates arrays must be supplied if semicoarsen orientation not given.");
      x = fineCoords->getDataNonConst(0);
      y = fineCoords->getDataNonConst(1);
      z = fineCoords->getDataNonConst(2);
      xptr = x.getRawPtr();
      yptr = y.getRawPtr();
      zptr = z.getRawPtr();
    }
    Teuchos::ArrayRCP<LO>     TLayerId   = Teuchos::arcp<LO>(Nnodes+1);  LayerId   = TLayerId.getRawPtr();
    Teuchos::ArrayRCP<LO>     TVertLineId= Teuchos::arcp<LO>(Nnodes+1);  VertLineId= TVertLineId.getRawPtr();

    NumZDir = ML_compute_line_info(LayerId, VertLineId,Ndofs, BlkSize,
                                   Zorientation, NumZDir,xptr,yptr,zptr, *(rowMap->getComm()));

    RCP<const Map> theCoarseMap;
    RCP<Matrix>    P;

    Ncoarse = MakeSemiCoarsenP(Nnodes,NumZDir,CoarsenRate,LayerId,VertLineId,
                               BlkSize, A, P, theCoarseMap);

    Teuchos::ArrayRCP<LO> coarseSemiInfo = Teuchos::arcp<LO>(3);
    coarseSemiInfo[NUM_ZPTS] = NumZDir*Ncoarse/Ndofs;
    coarseSemiInfo[ORIENTATION] = VERTICAL;
    coarseLevel.Set("SemiCoarsenInfo", coarseSemiInfo, MueLu::NoFactory::get());  // rst: why do I need NoFactory?
    Set(coarseLevel, "CoarseMap", P->getDomainMap());
    Set(coarseLevel, "P", P);

    // rst: null space might get scaled here ... do we care. We could just inject at the cpoints, but I don't
    //  feel that this is needed.

    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(), fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);

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

    if (PtsPerLine == 1) { printf("cannot coarsen further\n"); return -1; }
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
  LocalOrdinal SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MakeSemiCoarsenP(LocalOrdinal const Ntotal, LocalOrdinal const nz, LocalOrdinal const CoarsenRate, LocalOrdinal const LayerId[],
                                                                                                             LocalOrdinal const VertLineId[], LocalOrdinal const DofsPerNode, RCP<Matrix> & Amat, RCP<Matrix>& P, RCP<const Map>& coarseMap) const {

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
    int    BlkRow, dof_i, dof_j, node_k, *Sub2FullMap=NULL, RowLeng;
    int    i, j, iii, col, count, index, loc, PtRow, PtCol;
    SC *BandSol=NULL, *BandMat=NULL, TheSum;
    int    *IPIV=NULL, KL, KU, KLU, N, NRHS, LDAB,INFO;
    int    *Pcols;
    size_t *Pptr;
    SC *Pvals;
    int    MaxStencilSize, MaxNnzPerRow;
    int    *LayDiff=NULL;
    int    CurRow, LastGuy = -1, NewPtr;
    int    Ndofs;
    int    Nghost;
    int    *Layerdofs = NULL, *Col2Dof = NULL;

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

    RCP<Vector> localdtemp = VectorFactory::Build(Amat->getDomainMap());
    RCP<Vector> dtemp      = VectorFactory::Build(Amat->getColMap());
    ArrayRCP<SC> valptr= localdtemp->getDataNonConst(0);

    for (i = 0; i < Ntotal*DofsPerNode; i++)
      valptr[i]= (SC)LayerId[i/DofsPerNode];

    RCP< const Import> importer;
    importer = Amat->getCrsGraph()->getImporter();
    if (importer == Teuchos::null) {
      importer = ImportFactory::Build(Amat->getDomainMap(), Amat->getColMap());
    }
    dtemp->doImport(*localdtemp, *(importer), Xpetra::INSERT);

    valptr= dtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Layerdofs[i]= Teuchos::as<LO>( valptr[i] );
    valptr= localdtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode;        i++) valptr[i]= (SC) (i%DofsPerNode);
    dtemp->doImport(*localdtemp, *(importer), Xpetra::INSERT);
    valptr= dtemp->getDataNonConst(0);
    for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Col2Dof[i]= Teuchos::as<LO>( valptr[i] );

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
    NRHS = DofsPerNode;
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
    coarseMap = MapFactory::createUniformContigMapWithNode(rowMap->lib(),NCLayers*NVertLines*DofsPerNode,(rowMap->getComm()), rowMap->getNode());
    P       = rcp(new CrsMatrixWrap(rowMap, coarseMap , 0, Xpetra::StaticProfile));


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
            for (dof_i=0; dof_i < DofsPerNode; dof_i++) {

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
            for (dof_i = 0; dof_i < DofsPerNode; dof_i++) {
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
          for (dof_i=0; dof_i < DofsPerNode; dof_i++) {
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_compute_line_info(LocalOrdinal LayerId[], LocalOrdinal VertLineId[], LocalOrdinal Ndof, LocalOrdinal DofsPerNode, LocalOrdinal MeshNumbering, LocalOrdinal NumNodesPerVertLine, Scalar *xvals, Scalar *yvals, Scalar *zvals, const Teuchos::Comm<int>& comm) const {

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
        printf("Warning: did not assign %d to a vertical line?????\n",i);
      }
      if (LayerId[i] == -1) {
        printf("Warning: did not assign %d to a Layer?????\n",i);
      }
    }
    maxAll(&comm, NumNodesPerVertLine, i);
    if (NumNodesPerVertLine == -1)  NumNodesPerVertLine = i;

    TEUCHOS_TEST_FOR_EXCEPTION(NumNodesPerVertLine != i,Exceptions::RuntimeError, "Different processors have different z direction line lengths?\n");
    return NumNodesPerVertLine;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ML_az_dsort2(Scalar dlist[], LocalOrdinal N, LocalOrdinal list2[]) const {
    LO l, r, j, i, flag;
    LO RR2;
    SC       dRR, dK;

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
              if (dlist[j] > dlist[j - 1]) j = j + 1;

            if (dlist[j - 1] > dK) {
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
              if (dlist[j] > dlist[j - 1]) j = j + 1;
            if (dlist[j - 1] > dK) {
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

#define MUELU_SEMICOARSENPFACTORY_SHORT
#endif // MUELU_SEMICOARSENPFACTORY_DEF_HPP
