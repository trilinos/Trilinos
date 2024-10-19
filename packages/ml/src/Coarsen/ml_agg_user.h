/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_AGG_USER_H
#define ML_AGG_USER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int ML_SetUserLabel(char *user());

int ML_SetUserNumAggr(int (user)(ML_Operator*));

int ML_SetUserPartitions(int (user)(ML_Operator* Amat, char* bdry_nodes,
                                    double epsilon,
                                    double* x,double* y,double* z,
                                    int* partitions, int* LocalNonzeros));

extern int ML_Aggregate_CoarsenUser( ML_Aggregate *ml_ag,
                                    ML_Operator *Amatrix,
                                    ML_Operator **Pmatrix, ML_Comm *comm);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_AGG_USER_H */
