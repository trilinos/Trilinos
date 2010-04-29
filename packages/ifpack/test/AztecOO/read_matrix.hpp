#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <fstream>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

Epetra_CrsMatrix*
read_matrix_mm(const std::string& mm_file,
               const Epetra_Comm& comm);

Epetra_Vector*
read_vector_mm(const std::string& mm_file,
               const Epetra_Comm& comm);

void read_matrix_hb(const std::string& hb_file,
                    const Epetra_Comm& Comm,
                    Epetra_CrsMatrix*& A,
                    Epetra_Vector*& b);
#endif
