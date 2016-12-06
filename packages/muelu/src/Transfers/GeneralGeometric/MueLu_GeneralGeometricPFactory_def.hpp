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
#ifndef MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP
#define MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>


#include <Teuchos_LAPACK.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_GeneralGeometricPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // Coarsen can come in two forms, either a single char that will be interpreted as an integer which is used as the coarsening
    // rate in every spatial dimentions, or it can be a longer string that will then be interpreted as an array of integers.
    // By default coarsen is set as "2", hence a coarsening rate of 2 in every spatial dimension is the default setting!
    validParamList->set< std::string>            ("Coarsen",                   "2", "Coarsening rate in all spatial dimensions");
    validParamList->set< RCP<const FactoryBase> >("A",               Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",       Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",     Teuchos::null, "Generating factory for coorindates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "Coordinates");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    // obtain general variables
    using xdMV = Xpetra::MultiVector<double,LO,GO,NO>;
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    RCP<xdMV> fineCoords = Get< RCP<xdMV> >(fineLevel, "Coordinates");
    RCP<xdMV> coarseCoords;

    // get user-provided coarsening rate parameter (constant over all levels)
    const ParameterList& pL = GetParameterList();

    // collect general input data
    LO blkSize            = A->GetFixedBlockSize();
    RCP<const Map> rowMap = A->getRowMap();
    LO ndofs              = rowMap->getNodeNumElements();
    LO nNodes             = ndofs/blkSize;
    LO ndm                = 0;
    LO nFpts              = 0;

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords==Teuchos::null, Exceptions::RuntimeError, "Coordinates cannot be accessed from fine level!");
    ndm = fineCoords->getNumVectors();        // ndm stores the number of spatial dimensions
    nFpts = fineCoords->getGlobalLength();    // nFpts stores the number of points on the fine level


    // Get the coarsening rate
    std::string coarsenRate = pL.get<std::string>("Coarsen");
    ArrayRCP<LO> coarseRate = arcp<LO>(3);
    LO *coarse_rate = coarseRate.getRawPtr();
    if(coarsenRate.length()==1) {
      LO crate=Teuchos::as<LO>(Teuchos::ValueTypeConversionTraits< int, std::string >::safeConvert(coarsenRate));
      for(LO i = 0; i < 3; ++i) {
	if(i < ndm) {
	  coarse_rate[i] = crate;
	} else {
	  coarse_rate[i] = 1;
	}
      }
    } else {
      Teuchos::Array<LO> crates = Teuchos::fromStringToArray<int>(coarsenRate);
      TEUCHOS_TEST_FOR_EXCEPTION(crates.size()<ndm, Exceptions::RuntimeError, "Coarsen as a vector must have at least as many components as there are spatial dimensions in the problem.");
      for(LO i = 0; i < 3; ++i) {
	if(i < ndm) {
	  coarse_rate[i] = crates[i];
	} else {
	  coarse_rate[i] = 1;
	}
      }
    }

    /* Here one element can represent either the degenerate case of one node or the more general case of two nodes.
       i.e. x---x is a 1D element with two nodes and x is a 1D element with one node. This helps generating a 3D space from tensorial products... */
    LO nCpts = 1;
    LO neh[3];
    LO neH[3];
    LO end_rate[3];
    LO nTerms = 1;
    for(LO i = 0; i < 3; ++i) {
      if(i<ndm) {
	neh[i] = ceil(std::pow(nFpts,1.0/ndm)) - 1;
	end_rate[i] = neh[i] - coarse_rate[i]*(neh[i]/coarse_rate[i]);

	// Check for end_rate==0, which really means that nothing special needs to be done
	if(end_rate[i]==0) {
	  neH[i]  = (neh[i]/coarse_rate[i]);
	  end_rate[i] = coarse_rate[i];
	} else {
	  neH[i] = (neh[i]/coarse_rate[i]) + 1;
	}
      } else {
	neh[i] = 0;
	end_rate[i] = 1;
	neH[i] = 0;
      }


      if(i < ndm) {
	// From all this data we can compute the number of coarse points NCpts and the number of coefficients in the prolongation matrix.
	// 1D global  o---x---x ... x---x---+--o
	// left "o" has coarse_rate term
	// each "x" has 2*coarse_rate-1 terms and there are (neH+1)-3 "x"
	// "+" has coarse_rate+end_rate-1 terms
	// right "o" has end_rate terms
	//
	// So the grand total is Nterms = coarse_rate + ((neH+1)-3)*(2*coarse_rate-1) + coarse_rate+end_rate-1 + end_rate
	//                              = (neH-1)*(2*coarse_rate-1) + 2*end_rate
	nCpts = nCpts*(neH[i] + 1);
	nTerms = nTerms*((neH[i] - 1)*(2*coarse_rate[i] - 1) + 2*end_rate[i]);
      }
    }
    // For systems of PDEs we assume that all dofs have the same P operator.
    nTerms=nTerms*blkSize;

    // All that is left to do is loop over NCpts and:
    //   - extract coarse points coordiate for coarseCoords
    //   - get coordinates for current stencil computation
    //   - compute current stencil
    //   - compute row and column indices for stencil entries
    LO cNodePerDir[3];
    LO fNodePerDir[3];
    for(LO i = 0; i < 3; ++i) {
      cNodePerDir[i] = neH[i] + 1;
      fNodePerDir[i] = neh[i] + 1;
    }
    RCP<const Map> theCoarseMap;
    RCP<Matrix>    P;
    LO Ncoarse = MakeGeneralGeometricP(ndm,fNodePerDir,cNodePerDir,coarse_rate,end_rate,
				       fineCoords, nTerms, blkSize, theCoarseMap, A, P,
				       coarseCoords);

    // set StridingInformation of P   <---- what is happening here?
    if (A->IsView("stridedMaps") == true)
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), theCoarseMap);
    else
    P->CreateView("stridedMaps", P->getRangeMap(), theCoarseMap);

    // store the transfer operator and node coordinates on coarse level
    Set(coarseLevel, "P", P);
    Set(coarseLevel, "coarseCoordinates", coarseCoords);

    // rst: null space might get scaled here ... do we care. We could just inject at the cpoints, but I don't
    //  feel that this is needed.
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(), fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MakeGeneralGeometricP(LocalOrdinal const ndm, LocalOrdinal const *nFpts, LocalOrdinal const *nCpts, LocalOrdinal const coarse_rate[3], LocalOrdinal const end_rate[3], const RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> >& fCoords, LocalOrdinal const nnzP, LocalOrdinal const dofsPerNode, RCP<const Map>& coarseMap, RCP<Matrix> & Amat, RCP<Matrix>& P, RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> >& cCoords) const {

    /*
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * On termination, return the number of local prolongator columns owned by
     * this processor.
     *
     * Input
     * =====
     *    nNodes       Number of fine level Blk Rows owned by this processor
     *    coarse_rate  Rate of coarsening in each spatial direction.
     *    end_rate     Rate of coarsening in each spatial direction for the last
     *                 nodes in the mesh where an adaptive coarsening rate is
     *                 required.
     *    nTerms       Number of nonzero entries in the prolongation matrix.
     *    dofsPerNode  Number of degrees-of-freedom per mesh node.
     *
     * Output
     * =====
     *    So far nothing...
     */

    typedef Teuchos::ScalarTraits<SC> STS;
    const SC     zero      = STS::zero();
    const SC     one       = STS::one();
    const LO     INVALID   = Teuchos::OrdinalTraits<LO>::invalid();

    // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    // Teuchos::FancyOStream& out = *fancy;

    LO nFineNodes = nFpts[0]*nFpts[1]*nFpts[2];
    LO numRows = dofsPerNode*nFineNodes;
    LO nCoarseNodes = nCpts[0]*nCpts[1]*nCpts[2];
    LO numGloCols = dofsPerNode*nCoarseNodes;
    LO numLocCols = dofsPerNode*nCoarseNodes; // This Needs to change for parallel implementation!
    LO maxEntriesPerRow = 0;
    if(ndm==1) {
      maxEntriesPerRow = 2;
    } else if(ndm==2) {
      maxEntriesPerRow = 4;
    } else if(ndm==3) {
      maxEntriesPerRow = 8;
    }

    RCP<const Map> rowMap = Amat->getRowMap();
    std::vector<size_t> stridingInfo_;
    stridingInfo_.push_back(dofsPerNode);
    coarseMap=Xpetra::StridedMapFactory<LO,GO,NO>::Build(rowMap->lib(), Teuchos::as<size_t>(numGloCols),
							 Teuchos::as<size_t>(numLocCols), 0,
							 stridingInfo_, rowMap->getComm());
    P = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

    ArrayRCP<size_t>  iaP;
    ArrayRCP<LO>      jaP;
    ArrayRCP<SC>     valP;

    PCrs->allocateAllValues(nnzP, iaP, jaP, valP);

    ArrayView<size_t> ia  = iaP();
    ArrayView<LO>     ja  = jaP();
    ArrayView<SC>     val = valP();

    LO nMaxStencil = 0; LO nStencil = 0;
    LO iH = 0; LO count=0; LO ifH = 0;
    LO indices[4][3];
    std::div_t tmp;
    ia[0] = 0;

    // Declaration and assignment of fineCoords which holds the coordinates of the fine nodes in 3D.
    // To do so we pull the nD coordinates from fCoords and pad the rest with zero vectors...
    // We also create an nD array that will store the coordinates of the coarse grid nodes.
    // That array will eventually be used during the construction of cCoords, the MultiVector of
    // coarse nodes coordinates that we store on the coarse level.
    ArrayRCP< ArrayRCP<double> > fineCoords(3);
    RCP<const Map> fCoordsdMap=fCoords->getMap();
    RCP<Xpetra::Vector<double, LO, GO, NO> > zeros = Xpetra::VectorFactory<double, LO, GO, NO>::Build(fCoordsdMap, true);

    // Build the MultiVector holding the coarse grid points coordinates.
    // Parallel: MapFactory::Build will need to use the local number of coarse nodes instead of nCoarseNodes eventually!
    RCP<const Map> cCoordsMap=MapFactory::Build (fCoords->getMap()->lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
						 Teuchos::as<size_t>(nCoarseNodes),fCoords->getMap()->getIndexBase(),fCoords->getMap()->getComm());
    cCoords = Xpetra::MultiVectorFactory<double,LO,GO,Node>::Build(cCoordsMap,Teuchos::as<size_t>(ndm));
    ArrayRCP<double> xCoarseNodes; ArrayRCP<double> yCoarseNodes; ArrayRCP<double> zCoarseNodes;

    if(ndm==1) {
      fineCoords[0] = fCoords->getDataNonConst(0);
      fineCoords[1] = zeros->getDataNonConst(0);
      fineCoords[2] = zeros->getDataNonConst(0);
      nMaxStencil = 2;

      xCoarseNodes = cCoords->getDataNonConst(0);
    } else if(ndm==2) {
      fineCoords[0] = fCoords->getDataNonConst(0);
      fineCoords[1] = fCoords->getDataNonConst(1);
      fineCoords[2] = zeros->getDataNonConst(0);
      nMaxStencil = 4;

      xCoarseNodes=cCoords->getDataNonConst(0);
      yCoarseNodes=cCoords->getDataNonConst(1);
    } else if(ndm==3) {
      fineCoords[0] = fCoords->getDataNonConst(0);
      fineCoords[1] = fCoords->getDataNonConst(1);
      fineCoords[2] = fCoords->getDataNonConst(2);
      nMaxStencil = 8;

      xCoarseNodes = cCoords->getDataNonConst(0);
      yCoarseNodes = cCoords->getDataNonConst(1);
      zCoarseNodes = cCoords->getDataNonConst(2);
    }

    // First things: find the node ordering assuming prism element
    LO ordering=0;  //0 -> xyz; 1 -> xzy; 2 -> yxz; 3 -> yzx; 4 -> zxy; 5 -> zyx;
    if((fineCoords[1][0]==fineCoords[1][1]) && (fineCoords[2][0]==fineCoords[2][1])) {
      if(ndm==1) {
	ordering = 0;
      } else if(fineCoords[2][nFpts[0]-1]==fineCoords[2][nFpts[0]]) {
	ordering = 0;
      } else {
	ordering = 1;
      }
    } else if((fineCoords[0][0]==fineCoords[0][1]) && (fineCoords[2][0]==fineCoords[2][1])) {
      if(ndm==1) {
	ordering = 2;
      } else if(fineCoords[2][nFpts[1]-1]==fineCoords[2][nFpts[1]]) {
	ordering = 2;
      } else {
	ordering = 3;
      }
    } else {
      if(ndm==1) {
	ordering = 4;
      } else if(fineCoords[1][nFpts[2]-1]==fineCoords[1][nFpts[2]]) {
	ordering = 4;
      } else {
	ordering = 5;
      }
    }

    // The following array stores column offsets that are used to fill P
    const LO col_offsets[8] = {0, 1, nCpts[0], nCpts[0]+1, nCpts[1]*nCpts[0], nCpts[1]*nCpts[0]+1, nCpts[1]*nCpts[0]+nCpts[0], nCpts[1]*nCpts[0]+nCpts[0]+1};

    //Loop over the fine nodes and compute their interpolation stencils
    for(LO i = 0; i < nFineNodes; ++i) {
      LO nzIndStencil[8]={0,0,0,0,0,0,0,0};
      // Get point indices on fine grid
      tmp=std::div(i,nFpts[0]*nFpts[1]);
      indices[0][2] = tmp.quot;
      tmp=std::div(tmp.rem,nFpts[0]);
      indices[0][1] = tmp.quot;
      indices[0][0] = tmp.rem;

      // Get indices of ref point on coarse grid and location "flags"
      tmp=std::div(indices[0][0],coarse_rate[0]);
      indices[1][0] = tmp.quot;
      indices[2][0] = tmp.rem;
      tmp=std::div(indices[0][1],coarse_rate[1]);
      indices[1][1] = tmp.quot;
      indices[2][1] = tmp.rem;
      tmp=std::div(indices[0][2],coarse_rate[2]);
      indices[1][2] = tmp.quot;
      indices[2][2] = tmp.rem;

      // Get indices of ref point on fine grid
      indices[3][0] = indices[1][0]*coarse_rate[0];
      indices[3][1] = indices[1][1]*coarse_rate[1];
      indices[3][2] = indices[1][2]*coarse_rate[2];

      if(indices[0][0]==nFpts[0]-1) {
      	indices[1][0] = nCpts[0]-1;
      	indices[2][0] = 0;
      	indices[3][0] = nFpts[0]-1;
      }
      if(indices[0][1]==nFpts[1]-1) {
      	indices[1][1] = nCpts[1]-1;
      	indices[2][1] = 0;
      	indices[3][1] = nFpts[1]-1;
      }
      if(indices[0][2]==nFpts[2]-1) {
      	indices[1][2] = nCpts[2]-1;
      	indices[2][2] = 0;
      	indices[3][2] = nFpts[2]-1;
      }

      // Get reference point global index on coarse and fine grids
      iH =0;
      ifH=0;
      if(indices[3][0]==nFpts[0]-1) {
	ifH+= indices[3][0]-end_rate[0];
	iH += indices[1][0]-1;
      } else {
	ifH+= indices[3][0];
	iH += indices[1][0];
      }
      if(ndm > 1) {
	if(indices[3][1]==nFpts[1]-1) {
	  ifH+= (indices[3][1]-end_rate[1])*nFpts[0];
	  iH += (indices[1][1]-1)*nCpts[0];
	} else {
	  ifH+= indices[3][1]*nFpts[0];
	  iH += indices[1][1]*nCpts[0];
	}
      }
      if(ndm > 2) {
	if(indices[3][2]==nFpts[2]-1) {
	  ifH+= (indices[3][2]-end_rate[2])*nFpts[0]*nFpts[1];
	  iH += (indices[1][2]-1)*nCpts[0]*nCpts[1];
	} else {
	  ifH+= indices[3][2]*nFpts[0]*nFpts[1];
	  iH += indices[1][2]*nCpts[0]*nCpts[1];
	}
      }

      // Assuming lexicographic indexing the coarse nodes forming a prism
      // around fine node "i" are selected and store them in connec.
      // Two tricky things to be careful about:
      //    - are we using coarse_rate or end_rate?
      //      --> check indices and set rate correctly
      //    - are we on the east, north or top face?
      //      --> if so fix ifH to make sure there is
      //          a node to the east, north and top of ifH
      LO rate[3];
      if(indices[1][0] >= nCpts[0]-end_rate[0]-1) {
	rate[0] = end_rate[0];
      } else {
	rate[0] = coarse_rate[0];
      }
      if(indices[1][1] >= nCpts[1]-end_rate[1]-1) {
	rate[1] = end_rate[1];
      } else {
	rate[1] = coarse_rate[1];
      }
      if(indices[1][2] >= nCpts[2]-end_rate[2]-1) {
	rate[2] = end_rate[2];
      } else {
	rate[2] = coarse_rate[2];
      }
      if(ndm < 3) rate[2] = 0;
      if(ndm < 2) rate[1] = 0;

      SC connec[9][3];
      for(LO dim = 0; dim < 3; ++dim){
	connec[0][dim] = fineCoords[dim][i];
	connec[1][dim] = fineCoords[dim][ifH];
	connec[2][dim] = fineCoords[dim][ifH + rate[0]];
	connec[3][dim] = fineCoords[dim][ifH + rate[1]*nFpts[0]];
	connec[4][dim] = fineCoords[dim][ifH + rate[1]*nFpts[0] + rate[0]];
	connec[5][dim] = fineCoords[dim][ifH + rate[2]*nFpts[1]*nFpts[0]];
	connec[6][dim] = fineCoords[dim][ifH + rate[2]*nFpts[1]*nFpts[0] + rate[0]];
	connec[7][dim] = fineCoords[dim][ifH + rate[2]*nFpts[1]*nFpts[0] + rate[1]*nFpts[0]];
	connec[8][dim] = fineCoords[dim][ifH + rate[2]*nFpts[1]*nFpts[0] + rate[1]*nFpts[0] + rate[0]];
      }

      // Compute the actual geometric interpolation stencil
      SC stencil[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      ComputeStencil(ndm,ordering,connec,stencil);

      // Finally check whether the fine node is on a coarse: node, edge or face
      // and select accordingly the non-zero values from the stencil and the
      // corresponding column indices
      if(indices[2][0]==0 && indices[2][1]==0 && indices[2][2]==0) {
	nStencil = 1;
	// Fine node is on a coarse node, the local_node variable
	// tracks on which coarse node the fine node is located.
	LO cInd=indices[1][2]*nCpts[1]*nCpts[0]+indices[1][1]*nCpts[0]+indices[1][0];
	if(ndm==1) {
	  xCoarseNodes[cInd] = fineCoords[0][i];
	} else if(ndm==2) {
	  xCoarseNodes[cInd] = fineCoords[0][i];
	  yCoarseNodes[cInd] = fineCoords[1][i];
	} else if(ndm==3) {
	  xCoarseNodes[cInd] = fineCoords[0][i];
	  yCoarseNodes[cInd] = fineCoords[1][i];
	  zCoarseNodes[cInd] = fineCoords[2][i];
	}

	LO local_node=0; LO col_offset=0;
	if(indices[0][0]==nFpts[0]-1) {
	  nzIndStencil[0]+= 1;
	}
	if((indices[0][1]==nFpts[1]-1) && (ndm > 1)) {
	  nzIndStencil[0]+= 2;
	}
	if((indices[0][2]==nFpts[2]-1) && (ndm > 2)) {
	  nzIndStencil[0]+= 4;
	}
      } else if(indices[2][0]==0 && indices[2][1]==0) {
	// Fine node on z edge
	nStencil = 2;
	if((ndm > 2) && (indices[0][0]==nFpts[0]-1) && (indices[0][1]==nFpts[1]-1)) {
	  // Fine node on North-East edge
	  nzIndStencil[0] = 3;
	  nzIndStencil[1]  =7;
	} else if((ndm > 2) && (indices[0][1]==nFpts[1]-1)) {
	  // Fine node on North-West edge
	  nzIndStencil[0] = 2;
	  nzIndStencil[1] = 6;
	} else if((ndm > 1) && (indices[0][0]==nFpts[0]-1)) {
	  // Fine node on South-East edge
	  nzIndStencil[0] = 1;
	  nzIndStencil[1] = 5;
	} else {
	  // Fine node on South-West edge
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 4;
	}
      } else if(indices[2][0]==0 && indices[2][2]==0) {
	// Fine node on y edge
	nStencil = 2;
	if((ndm > 2) && (indices[0][0]==nFpts[0]-1) && (indices[0][2]==nFpts[2]-1)) {
	  // Fine node on Top-East edge
	  nzIndStencil[0] = 5;
	  nzIndStencil[1] = 7;
	} else if((ndm > 2) && (indices[0][2]==nFpts[2]-1)) {
	  // Fine node on Top-West edge
	  nzIndStencil[0] = 4;
	  nzIndStencil[1] = 6;
	} else if((ndm > 1) && (indices[0][0]==nFpts[0]-1)) {
	  // Fine node on Bottom-East edge
	  nzIndStencil[0] = 1;
	  nzIndStencil[1] = 3;
	} else {
	  // Fine node on Bottom-West edge
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 2;
	}
      } else if(indices[2][1]==0 && indices[2][2]==0) {
	// Fine node on x edge
	nStencil = 2;
	if((ndm > 2) && (indices[0][1]==nFpts[1]-1) && (indices[0][2]==nFpts[2]-1)) {
	  // Fine node on Top-North edge
	  nzIndStencil[0] = 6;
	  nzIndStencil[1] = 7;
	} else if((ndm > 2) && (indices[0][2]==nFpts[2]-1)) {
	  // Fine node on Top-South edge
	  nzIndStencil[0] = 4;
	  nzIndStencil[1] = 5;
	} else if((ndm > 1) && (indices[0][1]==nFpts[1]-1)) {
	  // Fine node on Bottom-North edge
	  nzIndStencil[0] = 2;
	  nzIndStencil[1] = 3;
	} else {
	  // Fine node on Bottom-South edge
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 1;
	}
      } else if(indices[2][0]==0) {
	// Fine node on yz face
	nStencil = 4;
	if((ndm > 2) && (indices[0][0]==nFpts[0]-1)) {
	  // Fine node on East face
	  nzIndStencil[0] = 1;
	  nzIndStencil[1] = 3;
	  nzIndStencil[2] = 5;
	  nzIndStencil[3] = 7;
	} else {
	  // Fine node on West face
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 2;
	  nzIndStencil[2] = 4;
	  nzIndStencil[3] = 6;
	}
      } else if(indices[2][1]==0) {
	// Fine node on xz face
	nStencil = 4;
	if((ndm > 2) && (indices[0][1]==nFpts[1]-1)) {
	  // Fine node on North face
	  nzIndStencil[0] = 2;
	  nzIndStencil[1] = 3;
	  nzIndStencil[2] = 6;
	  nzIndStencil[3] = 7;
	} else {
	  // Fine node on South face
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 1;
	  nzIndStencil[2] = 4;
	  nzIndStencil[3] = 5;
	}
      } else if(indices[2][2]==0) {
	// Fine node on xy face
	nStencil = 4;
	if((ndm > 2) && (indices[0][2]==nFpts[2]-1)) {
	  // Fine node on Top face
	  nzIndStencil[0] = 4;
	  nzIndStencil[1] = 5;
	  nzIndStencil[2] = 6;
	  nzIndStencil[3] = 7;
	} else {
	  // Fine node on Bottom face
	  nzIndStencil[0] = 0;
	  nzIndStencil[1] = 1;
	  nzIndStencil[2] = 2;
	  nzIndStencil[3] = 3;
	}
      } else {
	// Fine node in the interior of a hexahedron
	nStencil = 8;
	for(LO k = 0; k < 8; ++k) {
	  nzIndStencil[k] = k;
	}
      }

      // Here the values are filled in the Crs matrix arrays
      // This is basically the only place these variables are modified
      // Hopefully this makes handling system of PDEs easy!

      // Loop on dofsPerNode and process each row for the current Node
      for(LO j = 0; j < dofsPerNode; ++j) {
	ia[i*dofsPerNode + j + 1]=ia[i*dofsPerNode + j] + nStencil;
	for(LO k = 0; k < nStencil; ++k) {
	  ja[ia[i*dofsPerNode + j] + k] = (iH + col_offsets[nzIndStencil[k]])*dofsPerNode + j;
	  val[ia[i*dofsPerNode + j] + k] = stencil[nzIndStencil[k]];
	}
      }
    } // End loop over fine nodes

    // Set the values of the prolongation operators into the CrsMatrix P and call FillComplete
    PCrs->setAllValues(iaP, jaP, valP);
    PCrs->expertStaticFillComplete(coarseMap, Amat->getDomainMap());

    return 0;
  } // MakeGeneralGeometricP

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeStencil(const LO ndm, const LO ordering, const SC coord[9][3], SC stencil[8]) const {

    //                 ht    hn          7      8
    //                                   x------x
    //                  x  x            /|     /|
    //                  | /           5/ |   6/ |
    //                  |/            x------x  |
    //         hw  x----o----x  he    |  x---|--x
    //                 /|             | /3   | /4
    //                / |             |/     |/
    //  z  y         x  x             x------x
    //  | /                           1      2
    //  |/         hs    hb
    //  o---x
    //
    //  coords contains the coordinates of the node of interest indexed 0 and shown as "o" on the left plot
    //  and the coordinates of the hexahedron as indexed on the right plot.
    //  So far, assuming that the hexahedron is really a regular prism, the only coordinates
    //  needed are coord[0][], coord[1][] and coord[8][]
    //  hw=coord[0][0]-coord[1][0];  hs=coord[0][1]-coord[1][1];  hb=coord[0][2]-coord[1][2];
    //  he=coord[2][0]-coord[0][0];  hn=coord[2][1]-coord[0][1];  ht=coord[2][2]-coord[0][2];

    SC hw=0.0; SC he=0.0; SC hs=0.0; SC hn=0.0; SC hb=0.0; SC ht=0.0;
    hw = coord[0][0] - coord[1][0];
    he = coord[8][0] - coord[0][0];
    hs = coord[0][1] - coord[1][1];
    hn = coord[8][1] - coord[0][1];
    hb = coord[0][2] - coord[1][2];
    ht = coord[8][2] - coord[0][2];

    SC hx = 0.0; SC hy = 0.0; SC hz = 0.0;
    hx = hw + he;
    hy = hs + hn;
    hz = hb + ht;

    if(ndm==1) {
      hs = 1;
      hn = 1;
      hy = 1;
      hb = 1;
      ht = 1;
      hz = 1;
    } else if (ndm==2) {
      hb = 1;
      ht = 1;
      hz = 1;
    }

    LO permut[6][8] = {{0,1,2,3,4,5,6,7},
		       {0,1,4,5,2,3,6,7},
		       {0,2,1,3,4,6,5,7},
		       {0,4,1,5,2,6,3,7},
		       {0,2,4,6,1,3,5,7},
		       {0,4,2,6,1,5,3,7}};

    // The stencil needs reordering depending on the the node ordering: x-y-z or y-x-z, etc...
    stencil[permut[ordering][0]] = he/hx*hn/hy*ht/hz;
    stencil[permut[ordering][1]] = hw/hx*hn/hy*ht/hz;
    stencil[permut[ordering][2]] = he/hx*hs/hy*ht/hz;
    stencil[permut[ordering][3]] = hw/hx*hs/hy*ht/hz;
    stencil[permut[ordering][4]] = he/hx*hn/hy*hb/hz;
    stencil[permut[ordering][5]] = hw/hx*hn/hy*hb/hz;
    stencil[permut[ordering][6]] = he/hx*hs/hy*hb/hz;
    stencil[permut[ordering][7]] = hw/hx*hs/hy*hb/hz;
  } // End ComputeStencil
} //namespace MueLu

#define MUELU_GENERALGEOMETRICPFACTORY_SHORT
#endif // MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP
