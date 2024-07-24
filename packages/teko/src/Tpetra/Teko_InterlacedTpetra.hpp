// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_InterlacedTpetra_hpp__
#define __Teko_InterlacedTpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// include basic Tpetra information
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"

#include "Teko_ConfigDefs.hpp"

#include <vector>

namespace Teko {
namespace TpetraHelpers {
namespace Strided {

/** Build maps to make other conversions. This functions build maps assuming
 * that there are number of variables and you want to break them all up
 * into single blocks.
 *
 * \param[in] numGlobals The number of global unknowns in the original map
 * \param[in] numVars The number of sub maps to build
 * \param[in] comm Communicator to use in construction of the maps
 * \param[in,out] subMaps The maps for each block of unknowns. This vector
 *                        will be of length <code>numVars</code>. The integer
 *                        in the pair is the number of local variables being
 *                        represented. For this function it will always be 1.
 *                        The map itself will contain the global indices from
 *                        the original problem. This is useful for Export/Import
 *                        operations. The RCP for the map contains another map
 *                        associated with the string "contigMap" based on a
 *                        contiguous unknown numbering.
 *
 * \note <code> numGlobals % numVars == 0 </code>
 */
void buildSubMaps(GO numGlobals, int numVars, const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps);

/** Build maps to make other conversions. This functions build maps assuming
 * that there are number of variables and you want to break them all up
 * into blocks with sizes (in the number of varaibles) specified in <code>vars</code>
 *
 * \param[in] numGlobals The number of global unknowns in the original map
 * \param[in] vars Breakdown of each varaible to build.
 * \param[in] comm Communicator to use in construction of the maps
 * \param[in,out] subMaps The maps for each block of unknowns. This vector
 *                        will be of length <code>vars.size()</code>. The integer
 *                        in the pair is the number of local variables being
 *                        represented (this will be from <code>vars</code>).
 *                        The map itself will contain the global indices from
 *                        the original problem. This is useful for Export/Import
 *                        operations. The RCP for the map contains another map
 *                        associated with the string "contigMap" based on a
 *                        contiguous unknown numbering.
 *
 * \note <code> numGlobals % sum(vars) == 0 </code>
 */
void buildSubMaps(GO numGlobals, const std::vector<int>& vars, const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps);

/** Build maps to make other conversions. This functions build maps assuming
 * that there are number of variables and you want to break them all up
 * into blocks with sizes (in the number of varaibles) specified in <code>vars</code>
 *
 * \param[in] numGlobals The number of global unknowns in the original map
 * \param[in] numMyElements Number of nodes in this processor.
 * \param[in] minMyGID Minimum global ID on this processor.
 * \param[in] vars Breakdown of each varaible to build.
 * \param[in] comm Communicator to use in construction of the maps
 * \param[in,out] subMaps The maps for each block of unknowns. This vector
 *                        will be of length <code>vars.size()</code>. The integer
 *                        in the pair is the number of local variables being
 *                        represented (this will be from <code>vars</code>).
 *                        The map itself will contain the global indices from
 *                        the original problem. This is useful for Export/Import
 *                        operations. The RCP for the map contains another map
 *                        associated with the string "contigMap" based on a
 *                        contiguous unknown numbering.
 *
 * \note <code> numGlobals % sum(vars) == 0 </code>
 */
void buildSubMaps(GO numGlobals, LO numMyElements, GO minMyGID, const std::vector<int>& vars,
                  const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps);

/** Build maps to make other conversions. This functions build maps assuming
 * that there are number of variables and you want to break them all up
 * into blocks with sizes (in the number of varaibles) specified in <code>vars</code>
 *
 * \param[in] globalMap Parent map that describes the parallel layout of the new
 *                      sub maps.
 * \param[in] vars Breakdown of each varaible to build.
 * \param[in] comm Communicator to use in construction of the maps
 * \param[in,out] subMaps The maps for each block of unknowns. This vector
 *                        will be of length <code>vars.size()</code>. The integer
 *                        in the pair is the number of local variables being
 *                        represented (this will be from <code>vars</code>).
 *                        The map itself will contain the global indices from
 *                        the original problem. This is useful for Export/Import
 *                        operations. The RCP for the map contains another map
 *                        associated with the string "contigMap" based on a
 *                        contiguous unknown numbering.
 *
 * \note <code> numGlobals % sum(vars) == 0 </code>
 */
void buildSubMaps(const Tpetra::Map<LO, GO, NT>& globalMap, const std::vector<int>& vars,
                  const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps);

// build conversion import and export operators
void buildExportImport(
    const Tpetra::Map<LO, GO, NT>& baseMap,
    const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
    std::vector<Teuchos::RCP<Tpetra::Export<LO, GO, NT> > >& subExport,
    std::vector<Teuchos::RCP<Tpetra::Import<LO, GO, NT> > >& subImport);

// build a vector of subVectors
void buildSubVectors(
    const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
    std::vector<Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >& subVectors, LO count);

/** Associate a set of multi-vectors with a set of sub maps (this modifies the "extra data"
 * of the subVectors to contain info about the maps.
 */
void associateSubVectors(
    const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
    std::vector<Teuchos::RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >& subVectors);

// build a single subblock Epetra_CrsMatrix
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildSubBlock(
    int i, int j, const Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
    const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps);

// Rebuild a single subblock of a matrix
void rebuildSubBlock(
    int i, int j, const Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
    const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
    Tpetra::CrsMatrix<ST, LO, GO, NT>& mat);

// copy contents of many subvectors to a single vector
void many2one(Tpetra::MultiVector<ST, LO, GO, NT>& one,
              const std::vector<Teuchos::RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >& many,
              const std::vector<Teuchos::RCP<Tpetra::Export<LO, GO, NT> > >& subExport);

// copy contents of a single vector to many subvectors
void one2many(std::vector<Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >& many,
              const Tpetra::MultiVector<ST, LO, GO, NT>& single,
              const std::vector<Teuchos::RCP<Tpetra::Import<LO, GO, NT> > >& subImport);

}  // namespace Strided
}  // end namespace TpetraHelpers
}  // end namespace Teko

#endif
