#include "Teuchos_RCP.hpp"
#include "Epetra_CrsMatrix.h"
using namespace Teuchos;
RCP<Epetra_CrsMatrix> NewMatNewMap(Epetra_CrsMatrix& In, 
					   int Diagonal,
					   int ReindexRowMap,
					   int ReindexColMap,
					   int RangeMapType,
					   int DomainMapType
					   ) ;
