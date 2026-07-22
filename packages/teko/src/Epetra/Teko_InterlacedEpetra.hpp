// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_InterlacedEpetra_hpp__
#define __Teko_InterlacedEpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// include basic Epetra information
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include <vector>

namespace Teko {
namespace Epetra {
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
void buildSubMaps(int numGlobals, int numVars, const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps);

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
void buildSubMaps(int numGlobals, const std::vector<int>& vars, const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps);

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
void buildSubMaps(int numGlobals, int numMyElements, int minMyGID, const std::vector<int>& vars,
                  const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps);

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
void buildSubMaps(const Epetra_Map& globalMap, const std::vector<int>& vars,
                  const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps);

// build conversion import and export operators
void buildExportImport(const Epetra_Map& baseMap,
                       const std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps,
                       std::vector<Teuchos::RCP<Epetra_Export> >& subExport,
                       std::vector<Teuchos::RCP<Epetra_Import> >& subImport);

// build a vector of subVectors
void buildSubVectors(const std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps,
                     std::vector<Teuchos::RCP<Epetra_MultiVector> >& subVectors, int count);

/** Associate a set of multi-vectors with a set of sub maps (this modifies the "extra data"
 * of the subVectors to contain info about the maps.
 */
void associateSubVectors(const std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps,
                         std::vector<Teuchos::RCP<const Epetra_MultiVector> >& subVectors);

// build a single subblock Epetra_CrsMatrix
Teuchos::RCP<Epetra_CrsMatrix> buildSubBlock(
    int i, int j, const Epetra_CrsMatrix& A,
    const std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps);

// Rebuild a single subblock of a matrix
void rebuildSubBlock(int i, int j, const Epetra_CrsMatrix& A,
                     const std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps,
                     Epetra_CrsMatrix& mat);

// copy contents of many subvectors to a single vector
void many2one(Epetra_MultiVector& one,
              const std::vector<Teuchos::RCP<const Epetra_MultiVector> >& many,
              const std::vector<Teuchos::RCP<Epetra_Export> >& subExport);

// copy contents of a single vector to many subvectors
void one2many(std::vector<Teuchos::RCP<Epetra_MultiVector> >& many,
              const Epetra_MultiVector& single,
              const std::vector<Teuchos::RCP<Epetra_Import> >& subImport);

}  // namespace Strided
}  // end namespace Epetra
}  // end namespace Teko

#endif
