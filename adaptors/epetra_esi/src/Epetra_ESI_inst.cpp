
#include "Epetra_Array.h"
#include "Epetra_Array.cpp"

#include "Epetra_ESI_Argv.h"
#include "Epetra_ESI_Argv.cpp"

#include "Epetra_ESI_utils.h"
#include "Epetra_ESI_utils.cpp"

#include "Epetra_ESI_Object.h"
#include "Epetra_ESI_Object.cpp"

#include "Epetra_ESI_IndexSpace.h"
#include "Epetra_ESI_IndexSpace.cpp"

#include "Epetra_ESI_Vector.h"
#include "Epetra_ESI_Vector.cpp"

#include "Epetra_ESI_CrsMatrix.h"
#include "Epetra_ESI_CrsMatrix.cpp"

#include "Epetra_ESI_ostream.h"
#include "Epetra_ESI_ostream.cpp"

#include "Epetra_ESI_Operator.h"
#include "Epetra_ESI_Operator.cpp"

#include "Aztec_ESI_Translator.h"
#include "Aztec_ESI_Translator.cpp"

#include "Aztec_ESI_Solver.h"
#include "Aztec_ESI_Solver.cpp"

#include "Trilinos_ESI_Broker.h"
#include "Trilinos_ESI_Broker.cpp"

template class Epetra_Array<double>;
template class Epetra_Array<int>;
template class Epetra_Array<void*>;
template class Epetra_Array<char*>;
template class Epetra_Array<const char*>;

template class epetra_esi::IndexSpace<int>;

template class epetra_esi::Vector<double,int>;

template class epetra_esi::CrsMatrix<double,int>;

template class epetra_esi::Operator<double,int>;

template class aztecoo_esi::Solver<double,int>;

