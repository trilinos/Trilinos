#ifndef PANZER_ASSEMBLY_ENGINE_INARGS_HPP
#define PANZER_ASSEMBLY_ENGINE_INARGS_HPP

#include "Teuchos_RCP.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_Map;

namespace panzer {

  class LinearObjContainer;

  class AssemblyEngineInArgs {
    public:

    AssemblyEngineInArgs(const Teuchos::RCP<panzer::LinearObjContainer> & ghostedContainer,
                         const Teuchos::RCP<panzer::LinearObjContainer> & container)
       : ghostedContainer_(ghostedContainer), container_(container)
    { }

    Teuchos::RCP<panzer::LinearObjContainer> ghostedContainer_;
    Teuchos::RCP<panzer::LinearObjContainer> container_;

    double alpha;
    double beta;
    double time;
    bool evaluate_transient_terms;

  };

}

#endif
