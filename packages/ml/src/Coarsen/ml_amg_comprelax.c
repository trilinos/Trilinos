/* ******************************************************************** */
/* newly added to generate coarse grid based on compatible relaxation   */
/* -------------------------------------------------------------------- */

int ML_Coarsen_Dynamic(ML *ml, int level, int inc_or_dec)
{
   int         i, length, *indices;
   double      *initsol, *rhs, *sol, *res;
   ML_Operator *Amat;
   ML_Smoother *smoother;
   int         (*fun)(void *, int, double *, int, double *);


   Amat   = &(ml->Amat[level]);
   length = Amat->outvec_leng;
   indices= (int *)    malloc( length * sizeof(int) );
   initsol= (double *) malloc( length * sizeof(double) );
   sol    = (double *) malloc( length * sizeof(double) );
   rhs    = (double *) malloc( length * sizeof(double) );
   res    = (double *) malloc( length * sizeof(double) );
   ML_Smoother_Create( &smoother, &(ml->SingleLevel[level]) );
   fun    = ML_Smoother_SGS;
   ML_Smoother_Set(smoother, ML_INTERNAL, NULL, fun, NULL, 10, 1.0, NULL);
   for ( i = 0; i < length; i++ ) rhs[i] = 0.0;
   for ( i = 0; i < length; i++ ) initsol[i] = sol[i] = 1.0;

   ML_Smoother_Apply(smoother, length, sol, length, rhs, ML_ZERO);
/*
   ML_Operator_Apply(Amat, length, sol, length, res);
   for ( i = 0; i < length; i++ ) res[i] = rhs[i] - res[i];
*/
   for ( i = 0; i < length; i++ ) indices[i] = i;
   for ( i = 0; i < length; i++ ) sol[i] = sol[i] / initsol[i];
   ML_split_dsort( sol, length, indices, 10 );

   for ( i = 0; i < 10; i++ )
   printf("chosen coarse points %d = %6d \n", i, indices[i]);
   return 0;
}


