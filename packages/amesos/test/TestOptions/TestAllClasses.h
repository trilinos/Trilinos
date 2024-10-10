#include "Epetra_CrsMatrix.h"
#include <vector>

int TestAllClasses( const std::vector<std::string> AmesosClasses,
                     int EpetraMatrixType,
                    const std::vector<bool> AmesosClassesInstalled,
                    Epetra_CrsMatrix *& Amat,
                    const bool transpose,
                    const bool verbose,
                    const bool symmetric,
                    const int Levels,
                    const double Rcond,
                    int Diagonal,
                    int ReindexRowMap,
                    int ReindexColMap,
                    int RangeMapType,
                    int DomainMapType,
                    bool distribute,
                    const char *filename,
                    double &maxrelerror,
                    double &maxrelresidual,
                    int &NumTests) ;


#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
