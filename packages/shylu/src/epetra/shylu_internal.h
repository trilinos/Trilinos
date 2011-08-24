/* Internal functions in shylu */
#ifndef SHYLU_INTERNAL_H
#define SHYLU_INTERNAL_H

int create_matrices
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config    // i/p: library configuration
);

int extract_matrices
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config,   // i/p: library configuration
    bool insertValues       // true if we need to insert the values into the
                            // matrices, false if fillComplete is already called
                            // and we need to replace the local values
);

#endif //SHYLU_INTERNAL_H
