#include "ifp_c_wrappers.h"
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

struct AZ_IFPACK_STRUCT {                  
/* Data passing structure. This user */
/* defined data structure is used to pass information through   */
/* Aztec and back into the user's matrix-vector product. */
  int nr, nc;
  void *precon, *bmat;
  double *input_vector;
  int user_precon;
  void *user_aux_ptr;
};

typedef struct AZ_IFPACK_STRUCT AZ_IFPACK;

void AZ_ifpack_prec_create(double *x, double *b,
			   int *options, double *params,
			   int *proc_config,
			   AZ_MATRIX *Amat, AZ_PRECOND **Prec);

void AZ_ifpack_iterate(double *x, double *b,
               int *options, double *params, 
               double *status, int *proc_config,
               AZ_MATRIX *Amat);

void AZ_ifpack_precon(double x[], int *, int *,
                     double *, AZ_MATRIX *Amat, AZ_PRECOND *prec);

void AZ_ifpack_prec_destroy(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat, AZ_PRECOND *Prec);

void az2ifp_blockmatrix (void **bmat, AZ_MATRIX *Amat);
