#ifndef __PB_InterlacedEpetra_hpp__
#define __PB_InterlacedEpetra_hpp__

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

namespace PB {
namespace Epetra {

// build maps to make other conversions
void buildSubMaps(int numGlobals,int numVars,const Epetra_Comm & comm,
                  std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps);

// build maps to make other conversions
void buildSubMaps(int numGlobals,const std::vector<int> & vars,const Epetra_Comm & comm,
                  std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps);

// build conversion import and export operators
void buildExportImport(const Epetra_Map & baseMap, 
                       const std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps,
                       std::vector<Teuchos::RCP<Epetra_Export> > & subExport,
                       std::vector<Teuchos::RCP<Epetra_Import> > & subImport);

// build a vector of subVectors
void buildSubVectors(const std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps,
                     std::vector<Teuchos::RCP<Epetra_MultiVector> > & subVectors,int count);

// build a single subblock Epetra_CrsMatrix
Teuchos::RCP<Epetra_CrsMatrix> buildSubBlock(int i,int j,const Epetra_CrsMatrix & A,
                                             const std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps);

// copy contents of many subvectors to a single vector
void many2one(Epetra_MultiVector & one, const std::vector<Teuchos::RCP<const Epetra_MultiVector> > & many,
              const std::vector<Teuchos::RCP<Epetra_Export> > & subExport);

// copy contents of a single vector to many subvectors
void one2many(std::vector<Teuchos::RCP<Epetra_MultiVector> > & many,const Epetra_MultiVector & single,
              const std::vector<Teuchos::RCP<Epetra_Import> > & subImport);

} // end namespace Epetra
} // end namespace PB

#endif
