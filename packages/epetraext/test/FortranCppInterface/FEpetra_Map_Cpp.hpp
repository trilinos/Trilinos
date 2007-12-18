#include "FEpetra_Map.h"
#include "Teuchos_RCP.hpp"


#ifndef FEPETRA_MAP_CPP_HPP
#define FEPETRA_MAP_CPP_HPP


class Epetra_Map;


namespace FEpetra {


using Teuchos::RCP;


const RCP<const Epetra_Map> getMap( MapID mapID );


} // namespace FEpetra


#endif // FEPETRA_MAP_CPP_HPP
