// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_EXTRACTSUBMATRICES_DECL_HPP
#define _FROSCH_EXTRACTSUBMATRICES_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_Timers.h>

#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>

#include <KokkosKernels_Utils.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map);

    // ----------------------------------------------------------- //
    // split ExtractLocalSubdomainMatrix into symbolic / compute
    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Symbolic(RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix,        // input  : globalMatrix, re-distributed with map
                                              RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix);  // output : local submatrix

    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Compute(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > subdomainMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > repeatedMatrix);

    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Compute(RCP<      Import<LO,GO,NO> >    scatter,
                                             RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > subdomainMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > repeatedMatrix);
    // ----------------------------------------------------------- //

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map,
                                                                SC value);

    template <class SC,class LO,class GO,class NO>
    int UpdateLocalSubdomainMatrix(RCP<Matrix<SC,LO,GO,NO> > globalMatrix,
                                   RCP<Map<LO,GO,NO> > &map,
                                   RCP<Matrix<SC,LO,GO,NO> > &localSubdomainMatrix);

    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrices(RCP<const Matrix<SC,LO,GO,NO> > k,
                         ArrayView<GO> indI,
                         RCP<const Matrix<SC,LO,GO,NO> > &kII,
                         RCP<const Matrix<SC,LO,GO,NO> > &kIJ,
                         RCP<const Matrix<SC,LO,GO,NO> > &kJI,
                         RCP<const Matrix<SC,LO,GO,NO> > &kJJ);

    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrix(RCP<const Matrix<SC,LO,GO,NO> > k,
                       ArrayView<GO> indI,
                       RCP<const Matrix<SC,LO,GO,NO> > &kII);

    template <class LO,class GO,class NO>
    int BuildSubgraph(RCP<const CrsGraph<LO,GO,NO> > k,
                      ArrayView<GO> indI,
                      RCP<const CrsGraph<LO,GO,NO> > &kII);
}

#endif
