//
// file Util.h
//
// Additional function prototypes -- these 2 functions
// are not included in Trilinos_Util.h
//

class Epetra_Comm;


void read_hb_matrix(char *data_file, int MyPID,
	      int *N_global, int *n_nonzeros, 
	      double **val, int **bindx);



void distrib_msr_matrix(const Epetra_Comm & Comm,
	      int *N_global, int *n_nonzeros, 
           int *N_update, int **update, 
	      double **val, int **bindx);
