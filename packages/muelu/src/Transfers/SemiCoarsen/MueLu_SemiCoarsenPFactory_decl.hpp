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
#ifndef MUELU_SEMICOARSENPFACTORY_DECL_HPP
#define MUELU_SEMICOARSENPFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SemiCoarsenPFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

#define  VERTICAL       1
#define  HORIZONTAL     2
#define  GRID_SUPPLIED -1
#define  NUM_ZPTS       0
#define  ORIENTATION    1

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class SemiCoarsenPFactory : public PFactory {
#undef MUELU_SEMICOARSENPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    SemiCoarsenPFactory() { }

    //! Destructor.
    virtual ~SemiCoarsenPFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build (Level& fineLevel, Level& coarseLevel) const;
    void BuildP(Level& fineLevel, Level& coarseLevel) const;

    //@}

  private:
    LO FindCpts(LO const PtsPerLine, LO const CoarsenRate, LO const Thin, LO **LayerCpts) const;
    LO MakeSemiCoarsenP(LO const Ntotal, LO const nz, LO const CoarsenRate, LO const LayerId[],
                                  LO const VertLineId[], LO const DofsPerNode, RCP<Matrix>& Amat,
                                  RCP<Matrix>& P, RCP<const Map>& coarseMap) const;


    LO ML_compute_line_info(LO LayerId[], LO VertLineId[],
                                      LO Ndof, LO DofsPerNode,
                                      LO MeshNumbering, LO NumNodesPerVertLine,
                                      SC *xvals, SC *yvals, SC *zvals,
                                      const Teuchos::Comm<int>& comm ) const ;

    void ML_az_dsort2(SC dlist[], LO N, LO list2[]) const;

  }; //class SemiCoarsenPFactory

} //namespace MueLu

#define MUELU_SEMICOARSENPFACTORY_SHORT
#endif // MUELU_SEMICOARSENPFACTORY_DECL_HPP
