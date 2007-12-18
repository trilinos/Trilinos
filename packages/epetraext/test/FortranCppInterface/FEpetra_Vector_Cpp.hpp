#include "FEpetra_Vector.h"
#include "Teuchos_RCP.hpp"


#ifndef FEPETRA_VECTOR_CPP_HPP
#define FEPETRA_VECTOR_CPP_HPP


class Epetra_Vector;


namespace FEpetra {


using Teuchos::RCP;


const RCP<Epetra_Vector> getVector( VectorID vectorID );


} // namespace FEpetra


#endif // FEPETRA_VECTOR_CPP_HPP
