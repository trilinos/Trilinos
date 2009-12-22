#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <fstream>

#include "Epetra_CrsMatrix.h"

Epetra_CrsMatrix*
read_matrix_mm(const std::string& mm_file,
               const Epetra_Comm& comm);

#endif
