/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "mat_dh_private.h"
#include "Parser_dh.h"
#include "Hash_i_dh.h"
#include "Mat_dh.h"
#include "Mem_dh.h"
#include "Vec_dh.h"

#define IS_UPPER_TRI 97
#define IS_LOWER_TRI 98
#define IS_FULL      99
static int isTriangular(int m, int *rp, int *cval);

/* Instantiates Aout; allocates storage for rp, cval, and aval arrays;
   uses rowLengths[] and rowToBlock[] data to fill in rp[].
*/
static void mat_par_read_allocate_private(Mat_dh *Aout, int n, 
                                    int *rowLengths, int *rowToBlock);

/* Currently, divides (partitions)matrix by contiguous sections of rows.
   For future expansion: use metis.
*/
void mat_partition_private(Mat_dh A, int blocks, int *o2n_row, int *rowToBlock);


static void convert_triples_to_scr_private(int m, int nz, 
                                           int *I, int *J, double *A, 
                                           int *rp, int *cval, double *aval);

#if 0
#undef __FUNC__
#define __FUNC__ "mat_dh_print_graph_private"
void mat_dh_print_graph_private(int m, int beg_row, int *rp, int *cval, 
                    double *aval, int *n2o, int *o2n, Hash_i_dh hash, FILE* fp)
{
  START_FUNC_DH
  int i, j, row, col;
  double val;
  bool private_n2o = false;
  bool private_hash = false;

  if (n2o == NULL) {
    private_n2o = true;
    create_nat_ordering_private(m, &n2o); CHECK_V_ERROR;
    create_nat_ordering_private(m, &o2n); CHECK_V_ERROR;
  }
 
  if (hash == NULL) {
    private_hash = true;
    Hash_i_dhCreate(&hash, -1); CHECK_V_ERROR;
  }

  for (i=0; i<m; ++i) {
    row = n2o[i];
    for (j=rp[row]; j<rp[row+1]; ++j) {
      col = cval[j];
      if (col < beg_row || col >= beg_row+m) {
        int tmp = col;

        /* nonlocal column: get permutation from hash table */
        tmp = Hash_i_dhLookup(hash, col); CHECK_V_ERROR;
        if (tmp == -1) { 
          sprintf(msgBuf_dh, "beg_row= %i  m= %i; nonlocal column= %i not in hash table",
                                beg_row, m, col); 
          SET_V_ERROR(msgBuf_dh);
        } else {
          col = tmp;
        }
      } else {
        col = o2n[col];
      }

      if (aval == NULL) { 
        val = _MATLAB_ZERO_;
      } else {
        val = aval[j];
      }
      fprintf(fp, "%i %i %g\n", 1+row+beg_row, 1+col, val);
    }
  }

  if (private_n2o) {
    destroy_nat_ordering_private(n2o); CHECK_V_ERROR;
    destroy_nat_ordering_private(o2n); CHECK_V_ERROR;
  }

  if (private_hash) {
    Hash_i_dhDestroy(hash); CHECK_V_ERROR;
  }
  END_FUNC_DH
}

#endif


/* currently only for unpermuted */
#undef __FUNC__
#define __FUNC__ "mat_dh_print_graph_private"
void mat_dh_print_graph_private(int m, int beg_row, int *rp, int *cval, 
                    double *aval, int *n2o, int *o2n, Hash_i_dh hash, FILE* fp)
{
  START_FUNC_DH
  int i, j, row, col;
  bool private_n2o = false;
  bool private_hash = false;
  int *work = NULL;

  work = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;

  if (n2o == NULL) {
    private_n2o = true;
    create_nat_ordering_private(m, &n2o); CHECK_V_ERROR;
    create_nat_ordering_private(m, &o2n); CHECK_V_ERROR;
  }
 
  if (hash == NULL) {
    private_hash = true;
    Hash_i_dhCreate(&hash, -1); CHECK_V_ERROR;
  }

  for (i=0; i<m; ++i) {
    for (j=0; j<m; ++j) work[j] = 0;
    row = n2o[i];
    for (j=rp[row]; j<rp[row+1]; ++j) {
      col = cval[j];

      /* local column */
      if (col >= beg_row || col < beg_row+m) {
        col = o2n[col];
      } 

      /* nonlocal column: get permutation from hash table */
      else {
        int tmp = col;

        tmp = Hash_i_dhLookup(hash, col); CHECK_V_ERROR;
        if (tmp == -1) { 
          sprintf(msgBuf_dh, "beg_row= %i  m= %i; nonlocal column= %i not in hash table",
                                beg_row, m, col); 
          SET_V_ERROR(msgBuf_dh);
        } else {
          col = tmp;
        }
      } 

      work[col] = 1;
    }

    for (j=0; j<m; ++j) {
      if (work[j]) {
        fprintf(fp, " x ");
      } else {
        fprintf(fp, "   ");
      }
    }
    fprintf(fp, "\n");
  }

  if (private_n2o) {
    destroy_nat_ordering_private(n2o); CHECK_V_ERROR;
    destroy_nat_ordering_private(o2n); CHECK_V_ERROR;
  }

  if (private_hash) {
    Hash_i_dhDestroy(hash); CHECK_V_ERROR;
  }

  if (work != NULL) { FREE_DH(work); CHECK_V_ERROR; }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "create_nat_ordering_private"
void create_nat_ordering_private(int m, int **p)
{
  START_FUNC_DH
  int *tmp, i;

  tmp = *p = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) {
    tmp[i] = i;
  }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "destroy_nat_ordering_private"
void destroy_nat_ordering_private(int *p)
{
  START_FUNC_DH
  FREE_DH(p); CHECK_V_ERROR;
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "invert_perm"
void invert_perm(int m, int *pIN, int *pOUT)
{
  START_FUNC_DH
  int i;

  for (i=0; i<m; ++i) pOUT[pIN[i]] = i;
  END_FUNC_DH
}



/* only implemented for a single cpu! */
#undef __FUNC__
#define __FUNC__ "mat_dh_print_csr_private"
void mat_dh_print_csr_private(int m, int *rp, int *cval, double *aval, FILE* fp)
{
  START_FUNC_DH
  int i, nz = rp[m];

  /* print header line */
  fprintf(fp, "%i %i\n", m, rp[m]);

  /* print rp[] */
  for (i=0; i<=m; ++i) fprintf(fp, "%i ", rp[i]);
  fprintf(fp, "\n");

  /* print cval[] */
  for (i=0; i<nz; ++i) fprintf(fp, "%i ", cval[i]);
  fprintf(fp, "\n");

  /* print aval[] */
  for (i=0; i<nz; ++i) fprintf(fp, "%1.19e ", aval[i]);
  fprintf(fp, "\n");

  END_FUNC_DH
}


/* only implemented for a single cpu! */
#undef __FUNC__
#define __FUNC__ "mat_dh_read_csr_private"
void mat_dh_read_csr_private(int *mOUT, int **rpOUT, int **cvalOUT, 
                                            double **avalOUT, FILE* fp)
{
  START_FUNC_DH
  int i, m, nz, items;
  int *rp, *cval;
  double *aval;

  /* read header line */
  items = fscanf(fp,"%d %d",&m, &nz);
  if (items != 2) {
    SET_V_ERROR("failed to read header");
  } else {
    printf("mat_dh_read_csr_private:: m= %i  nz= %i\n", m, nz);
  }

  *mOUT = m;
  rp = *rpOUT = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  cval = *cvalOUT = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  aval = *avalOUT = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;

  /* read rp[] block */
  for (i=0; i<=m; ++i) {
    items = fscanf(fp,"%d", &(rp[i]));
    if (items != 1) {
      sprintf(msgBuf_dh, "failed item %i of %i in rp block", i, m+1);
      SET_V_ERROR(msgBuf_dh);
    }
  }

  /* read cval[] block */
  for (i=0; i<nz; ++i) {
    items = fscanf(fp,"%d", &(cval[i]));
    if (items != 1) {
      sprintf(msgBuf_dh, "failed item %i of %i in cval block", i, m+1);
      SET_V_ERROR(msgBuf_dh);
    }
  }

  /* read aval[] block */
  for (i=0; i<nz; ++i) {
    items = fscanf(fp,"%lg", &(aval[i]));
    if (items != 1) {
      sprintf(msgBuf_dh, "failed item %i of %i in aval block", i, m+1);
      SET_V_ERROR(msgBuf_dh);
    }
  }
  END_FUNC_DH
}

/*============================================*/
#define MAX_JUNK 200

#undef __FUNC__
#define __FUNC__ "mat_dh_read_triples_private"
void mat_dh_read_triples_private(int ignore, int *mOUT, int **rpOUT, 
                                   int **cvalOUT, double **avalOUT, FILE* fp)
{
  START_FUNC_DH
  int m, n, nz, items, i, j;
  int idx = 0;
  int *cval, *rp, *I, *J;
  double *aval, *A, v;
  char junk[MAX_JUNK];
  fpos_t fpos;

  /* skip over header */
  if (ignore && myid_dh == 0) {
    printf("mat_dh_read_triples_private:: ignoring following header lines:\n");
    printf("--------------------------------------------------------------\n");
    for (i=0; i<ignore; ++i) {
      fgets(junk, MAX_JUNK, fp);
      printf("%s", junk);
    }
    printf("--------------------------------------------------------------\n");
    if (fgetpos(fp, &fpos)) SET_V_ERROR("fgetpos failed!");
    printf("\nmat_dh_read_triples_private::1st two non-ignored lines:\n");
    printf("--------------------------------------------------------------\n");
    for (i=0; i<2; ++i) {
      fgets(junk, MAX_JUNK, fp);
      printf("%s", junk);
    }
    printf("--------------------------------------------------------------\n");
    if (fsetpos(fp, &fpos)) SET_V_ERROR("fsetpos failed!");
  }


if (feof(fp)) printf("trouble!");

  /* determine matrix dimensions */
  m=n=nz=0;
  while (!feof(fp)) {
    items = fscanf(fp,"%d %d %lg",&i,&j,&v);
    if (items != 3) {
      break;
    }
    ++nz;
    if (i > m) m = i;
    if (j > n) n = j;
  }

  if (myid_dh == 0) {
    printf("mat_dh_read_triples_private: m= %i  nz= %i\n", m, nz);
  }


  /* reset file, and skip over header again */
  rewind(fp);
  for (i=0; i<ignore; ++i) {
    fgets(junk, MAX_JUNK, fp);
  }

  /* error check for squareness */
  if (m != n) {
    sprintf(msgBuf_dh, "matrix is not square; row= %i, cols= %i", m, n);
    SET_V_ERROR(msgBuf_dh);
  }

  *mOUT = m;

  /* allocate storage */
  rp = *rpOUT = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  cval = *cvalOUT = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  aval = *avalOUT = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;

  I = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  J = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  A = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;

  /* read <row, col, value> triples into arrays */
  while (!feof(fp)) {
    items = fscanf(fp,"%d %d %lg",&i,&j,&v);
    if (items < 3) break;
    j--;
    i--;
    I[idx] = i;
    J[idx] = j;
    A[idx] = v;
    ++idx;
  }

  /* convert from triples to sparse-compressed-row storage */
  convert_triples_to_scr_private(m, nz, I, J, A, rp, cval, aval); CHECK_V_ERROR;

  /* if matrix is triangular */
  { int type;
    type = isTriangular(m, rp, cval); CHECK_V_ERROR;
    if (type == IS_UPPER_TRI) {
      printf("CAUTION: matrix is upper triangular; converting to full\n");
    } else if (type == IS_LOWER_TRI) {
      printf("CAUTION: matrix is lower triangular; converting to full\n");
    }

    if (type == IS_UPPER_TRI || type == IS_LOWER_TRI) {
      make_full_private(m, &rp, &cval, &aval); CHECK_V_ERROR;
    }
  }

  *rpOUT = rp;
  *cvalOUT = cval;
  *avalOUT = aval;

  FREE_DH(I); CHECK_V_ERROR;
  FREE_DH(J); CHECK_V_ERROR;
  FREE_DH(A); CHECK_V_ERROR;

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "convert_triples_to_scr_private"
void convert_triples_to_scr_private(int m, int nz, int *I, int *J, double *A, 
                                      int *rp, int *cval, double *aval)
{
  START_FUNC_DH
  int i;
  int *rowCounts;

  rowCounts = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) rowCounts[i] =   0;

  /* count number of entries in each row */
  for (i=0; i<nz; ++i) {
    int row = I[i];
    rowCounts[row] += 1;
  }

  /* prefix-sum to form rp[] */
  rp[0] = 0;
  for (i=1; i<=m; ++i) {
    rp[i] = rp[i-1] + rowCounts[i-1];
  }
  memcpy(rowCounts, rp, (m+1)*sizeof(int));

  /* write SCR arrays */
  for (i=0; i<nz; ++i) {
    int row = I[i];
    int col = J[i];
    double val = A[i];
    int idx = rowCounts[row];
    rowCounts[row] += 1;

    cval[idx] = col;
    aval[idx] = val;
  }


  FREE_DH(rowCounts); CHECK_V_ERROR;
  END_FUNC_DH
}


/*======================================================================
 * utilities for use in drivers that read, write, convert, and/or
 * compare different file types
 *======================================================================*/

void fix_diags_private(Mat_dh A);
void insert_missing_diags_private(Mat_dh A);

#undef __FUNC__
#define __FUNC__ "readMat"
void readMat(Mat_dh *Aout, char *ft, char *fn, int ignore)
{
  START_FUNC_DH
  bool makeStructurallySymmetric;
  bool fixDiags;
  *Aout = NULL;

  makeStructurallySymmetric = 
      Parser_dhHasSwitch(parser_dh, "-makeSymmetric");
  fixDiags = 
      Parser_dhHasSwitch(parser_dh, "-fixDiags");

  if (fn == NULL) {
    SET_V_ERROR("passed NULL filename; can't open for reading!");
  }

  if (!strcmp(ft, "csr")) 
  {
    Mat_dhReadCSR(Aout, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "trip")) 
  {
    Mat_dhReadTriples(Aout, ignore, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "ebin"))
  {
    Mat_dhReadBIN(Aout, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "petsc")) {
    sprintf(msgBuf_dh, "must recompile Euclid using petsc mode!");
    SET_V_ERROR(msgBuf_dh);
  }

  else 
  {
    sprintf(msgBuf_dh, "unknown filetype: -ftin %s", ft);
    SET_V_ERROR(msgBuf_dh);
  }

  if (makeStructurallySymmetric) {
    printf("\npadding with zeros to make structurally symmetric\n");
    Mat_dhMakeStructurallySymmetric(*Aout); CHECK_V_ERROR;
  }

  if ( (*Aout)->m == 0) {
    SET_V_ERROR("row count = 0; something's wrong!");
  }

  if (fixDiags) {
    fix_diags_private(*Aout); CHECK_V_ERROR;
  }

  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "fix_diags_private"
void fix_diags_private(Mat_dh A)
{
  START_FUNC_DH
  int i, j, m = A->m, *rp = A->rp, *cval = A->cval;
  double *aval = A->aval;
  bool insertDiags = false;

  /* verify that all diagonals are present */
  for (i=0; i<m; ++i) {
    bool isMissing = true;
    for (j=rp[i]; j<rp[i+1]; ++j) {
      if (cval[j] == i) {
        isMissing = false;
        break;
      }
    }
    if (isMissing) {
      insertDiags = true;
      break;
    }
  }

  if (insertDiags) {
    insert_missing_diags_private(A); CHECK_V_ERROR;
    rp = A->rp; 
    cval = A->cval;
    aval = A->aval;
  }

  /* set value of all diags to largest absolute value in each row */
  for (i=0; i<m; ++i) {
    double sum = 0;
    for (j=rp[i]; j<rp[i+1]; ++j) {
      sum = MAX(sum, fabs(aval[j]));
    }
    for (j=rp[i]; j<rp[i+1]; ++j) {
      if (cval[j] == i) {
        aval[j] = sum;
        break;
      }
    }
  }

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "insert_missing_diags_private"
void insert_missing_diags_private(Mat_dh A)
{
  START_FUNC_DH
  int *RP = A->rp, *CVAL = A->cval, m = A->m;
  int *rp, *cval;
  double *AVAL = A->aval, *aval;
  int i, j, nz = RP[m]+m;
  int idx = 0;

  rp = A->rp = (int *)MALLOC_DH((1+m)*sizeof(int)); CHECK_V_ERROR;
  cval = A->cval = (int *)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  aval = A->aval = (double *)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;
  rp[0] = 0;

  for (i=0; i<m; ++i) {
    bool isMissing = true;
    for (j=RP[i]; j<RP[i+1]; ++j) {
      cval[idx] = CVAL[j];
      aval[idx] = AVAL[j];
      ++idx;
      if (CVAL[j] == i) isMissing = false;
    }
    if (isMissing) {
      cval[idx] = i;
      aval[idx] = 0.0;
      ++idx;
    }
    rp[i+1] = idx;
  }
  
  FREE_DH(RP); CHECK_V_ERROR;
  FREE_DH(CVAL); CHECK_V_ERROR;
  FREE_DH(AVAL); CHECK_V_ERROR;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "readVec"
void readVec(Vec_dh *bout, char *ft, char *fn, int ignore)
{
  START_FUNC_DH
  *bout = NULL;

  if (fn == NULL) {
    SET_V_ERROR("passed NULL filename; can't open for reading!");
  }

  if (!strcmp(ft, "csr")  ||  !strcmp(ft, "trip")) 
  {
    Vec_dhRead(bout, ignore, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "ebin"))
  {
    Vec_dhReadBIN(bout, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "petsc")) {
    sprintf(msgBuf_dh, "must recompile Euclid using petsc mode!");
    SET_V_ERROR(msgBuf_dh);
  }

  else 
  {
    sprintf(msgBuf_dh, "unknown filetype: -ftin %s", ft);
    SET_V_ERROR(msgBuf_dh);
  }
  
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "writeMat"
void writeMat(Mat_dh Ain, char *ft, char *fn)
{
  START_FUNC_DH
  if (fn == NULL) {
    SET_V_ERROR("passed NULL filename; can't open for writing!");
  }

  if (!strcmp(ft, "csr")) 
  {
    Mat_dhPrintCSR(Ain, NULL, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "trip")) 
  {
    Mat_dhPrintTriples(Ain, NULL, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "ebin"))
  {
    Mat_dhPrintBIN(Ain, NULL, fn); CHECK_V_ERROR;
  } 


  else if (!strcmp(ft, "petsc")) {
    sprintf(msgBuf_dh, "must recompile Euclid using petsc mode!");
    SET_V_ERROR(msgBuf_dh);
  }

  else 
  {
    sprintf(msgBuf_dh, "unknown filetype: -ftout %s", ft);
    SET_V_ERROR(msgBuf_dh);
  }

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "writeVec"
void writeVec(Vec_dh bin, char *ft, char *fn)
{
  START_FUNC_DH
  if (fn == NULL) {
    SET_V_ERROR("passed NULL filename; can't open for writing!");
  }

  if (!strcmp(ft, "csr")  ||  !strcmp(ft, "trip")) 
  {
    Vec_dhPrint(bin, NULL, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "ebin"))
  {
    Vec_dhPrintBIN(bin, NULL, fn); CHECK_V_ERROR;
  } 

  else if (!strcmp(ft, "petsc")) {
    sprintf(msgBuf_dh, "must recompile Euclid using petsc mode!");
    SET_V_ERROR(msgBuf_dh);
  }

  else 
  {
    sprintf(msgBuf_dh, "unknown filetype: -ftout %s", ft);
    SET_V_ERROR(msgBuf_dh);
  }

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "isTriangular"
int isTriangular(int m, int *rp, int *cval)
{
  START_FUNC_DH
  int row, j;
  int type;
  bool type_lower = false, type_upper = false;

  if (np_dh > 1) {
    SET_ERROR(-1, "only implemented for a single cpu");
  }

  for (row=0; row<m; ++row) {
    for (j=rp[row]; j<rp[row+1]; ++j) {
      int col = cval[j];
      if (col < row) type_lower = true;
      if (col > row) type_upper = true;
    }
    if (type_lower && type_upper) break;
  }

  if (type_lower && type_upper) {
    type = IS_FULL;
  } else if (type_lower) {
    type = IS_LOWER_TRI;
  } else {
    type = IS_UPPER_TRI;
  }
  END_FUNC_VAL(type)
}

/*-----------------------------------------------------------------------------------*/

static void mat_dh_transpose_reuse_private_private(
                              bool allocateMem, int m, 
                              int *rpIN, int *cvalIN, double *avalIN,
                              int **rpOUT, int **cvalOUT, double **avalOUT);


#undef __FUNC__
#define __FUNC__ "mat_dh_transpose_reuse_private"
void mat_dh_transpose_reuse_private(int m, 
                              int *rpIN, int *cvalIN, double *avalIN,
                              int *rpOUT, int *cvalOUT, double *avalOUT)
{
  START_FUNC_DH
  mat_dh_transpose_reuse_private_private(false, m, rpIN, cvalIN, avalIN,
                                       &rpOUT, &cvalOUT, &avalOUT); CHECK_V_ERROR;
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "mat_dh_transpose_private"
void mat_dh_transpose_private(int m, int *RP, int **rpOUT,
                              int *CVAL, int **cvalOUT,
                              double *AVAL, double **avalOUT)
{
  START_FUNC_DH
  mat_dh_transpose_reuse_private_private(true, m, RP, CVAL, AVAL, 
                                       rpOUT, cvalOUT, avalOUT); CHECK_V_ERROR;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "mat_dh_transpose_private_private"
void mat_dh_transpose_reuse_private_private(bool allocateMem, int m, 
                              int *RP, int *CVAL, double *AVAL,
                              int **rpOUT, int **cvalOUT, double **avalOUT)
{
  START_FUNC_DH
  int *rp, *cval, *tmp;
  int i, j, nz = RP[m];
  double *aval;

  if (allocateMem) {
    rp = *rpOUT = (int *)MALLOC_DH((1+m)*sizeof(int)); CHECK_V_ERROR;
    cval = *cvalOUT = (int *)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
    if (avalOUT != NULL) {
      aval = *avalOUT = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;
    }
  } else {
    rp = *rpOUT;
    cval = *cvalOUT;
    if (avalOUT != NULL) aval = *avalOUT;
  }


  tmp = (int *)MALLOC_DH((1+m)*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<=m; ++i) tmp[i] = 0;

  for (i=0; i<m; ++i) {
    for (j=RP[i]; j<RP[i+1]; ++j) {
      int col = CVAL[j];
      tmp[col+1] += 1;
    }
  }
  for (i=1; i<=m; ++i) tmp[i] += tmp[i-1];
  memcpy(rp, tmp, (m+1)*sizeof(int));

  if (avalOUT != NULL) {
    for (i=0; i<m; ++i) {
      for (j=RP[i]; j<RP[i+1]; ++j) {
        int col = CVAL[j];
        int idx = tmp[col];
        cval[idx] = i;
        aval[idx] = AVAL[j];
        tmp[col] += 1;
      }
    }
  }

  else {
    for (i=0; i<m; ++i) {
      for (j=RP[i]; j<RP[i+1]; ++j) {
        int col = CVAL[j];
        int idx = tmp[col];
        cval[idx] = i;
        tmp[col] += 1;
      }
    }
  }

  FREE_DH(tmp); CHECK_V_ERROR;
  END_FUNC_DH
}

/*-----------------------------------------------------------------------------------*/

#undef __FUNC__
#define __FUNC__ "mat_find_owner"
int mat_find_owner(int *beg_rows, int *end_rows, int index)
{
  START_FUNC_DH
  int pe, owner = -1;

  for (pe=0; pe<np_dh; ++pe) {
    if (index >= beg_rows[pe] && index < end_rows[pe]) {
      owner = pe;
      break;
    }
  }

  if (owner == -1) {
    sprintf(msgBuf_dh, "failed to find owner for index= %i", index);
    SET_ERROR(-1, msgBuf_dh);
  }

  END_FUNC_VAL(owner)
}


#define AVAL_TAG 2
#define CVAL_TAG 3
void partition_and_distribute_private(Mat_dh A, Mat_dh *Bout);
void partition_and_distribute_metis_private(Mat_dh A, Mat_dh *Bout); 

#undef __FUNC__
#define __FUNC__ "readMat_par"
void readMat_par(Mat_dh *Aout, char *fileType, char *fileName, int ignore)
{
  START_FUNC_DH
  Mat_dh A = NULL;

  if (myid_dh == 0) {
    int tmp = np_dh;
    np_dh = 1;
    readMat(&A, fileType, fileName, ignore); CHECK_V_ERROR;
    np_dh = tmp;
  }

  if (np_dh == 1) {
    *Aout = A;
  } else {
    if (Parser_dhHasSwitch(parser_dh, "-metis")) {
      partition_and_distribute_metis_private(A, Aout); CHECK_V_ERROR;
    } else {
      partition_and_distribute_private(A, Aout); CHECK_V_ERROR;
    }
  }

  if (np_dh > 1 && A != NULL) {
    Mat_dhDestroy(A); CHECK_V_ERROR;
  }

 
  if (Parser_dhHasSwitch(parser_dh, "-printMAT")) {
    char xname[] = "A", *name = xname;
    Parser_dhReadString(parser_dh, "-printMat", &name);
    Mat_dhPrintTriples(*Aout, NULL, name); CHECK_V_ERROR;
    printf_dh("\n@@@ readMat_par: printed mat to %s\n\n", xname);
  }


  END_FUNC_DH
}

/* this is bad code! */
#undef __FUNC__
#define __FUNC__ "partition_and_distribute_metis_private"
void partition_and_distribute_metis_private(Mat_dh A, Mat_dh *Bout)
{
  START_FUNC_DH
  Mat_dh B = NULL;
  Mat_dh C = NULL;
  int i, m;
  int *rowLengths = NULL;
  int *o2n_row = NULL, *n2o_col = NULL, *rowToBlock = NULL;
  int *beg_row = NULL, *row_count = NULL;
  MPI_Request *send_req = NULL;
  MPI_Request *rcv_req = NULL;
  MPI_Status  *send_status = NULL;
  MPI_Status  *rcv_status = NULL;

  MPI_Barrier(comm_dh);
  printf_dh("@@@ partitioning with metis\n");

  /* broadcast number of rows to all processors */
  if (myid_dh == 0)  m = A->m;
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* broadcast number of nonzeros in each row to all processors */
  rowLengths = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  rowToBlock = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;

  if (myid_dh == 0) {
    int *tmp = A->rp;
    for (i=0; i<m; ++i) {
      rowLengths[i] = tmp[i+1] - tmp[i];
    }
  }
  MPI_Bcast(rowLengths, m, MPI_INT, 0, comm_dh);

  /* partition matrix */
  if (myid_dh == 0) {
    int idx = 0;
    int j;

    /* partition and permute matrix */
    Mat_dhPartition(A, np_dh, &beg_row, &row_count, &n2o_col, &o2n_row); ERRCHKA;
    Mat_dhPermute(A, n2o_col, &C); ERRCHKA;
 
    /* form rowToBlock array */
    for (i=0; i<np_dh; ++i) {
      for (j=beg_row[i]; j<beg_row[i]+row_count[i]; ++j) {
        rowToBlock[idx++] = i;
      }
    }
  }

  /* broadcast partitiioning information to all processors */
  MPI_Bcast(rowToBlock, m, MPI_INT, 0, comm_dh);

  /* allocate storage for local portion of matrix */
  mat_par_read_allocate_private(&B, m, rowLengths, rowToBlock); CHECK_V_ERROR;

  /* root sends each processor its portion of the matrix */
  if (myid_dh == 0) {
    int *cval = C->cval, *rp = C->rp;
    double *aval = C->aval;
    send_req = (MPI_Request*)MALLOC_DH(2*m*sizeof(MPI_Request)); CHECK_V_ERROR;
    send_status = (MPI_Status*)MALLOC_DH(2*m*sizeof(MPI_Status)); CHECK_V_ERROR;
    for (i=0; i<m; ++i) {
      int owner = rowToBlock[i];
      int count = rp[i+1]-rp[i];

      /* error check for empty row */
      if (! count) {
        sprintf(msgBuf_dh, "row %i of %i is empty!", i+1, m);
        SET_V_ERROR(msgBuf_dh);
      }

      MPI_Isend(cval+rp[i], count, MPI_INT, owner, CVAL_TAG, comm_dh, send_req+2*i);
      MPI_Isend(aval+rp[i], count, MPI_DOUBLE, owner, AVAL_TAG, comm_dh, send_req+2*i+1);
    }
  } 

  /* all processors receive their local rows */
  { int *cval = B->cval;
    int *rp = B->rp;
    double *aval = B->aval;
    m = B->m;

    rcv_req = (MPI_Request*)MALLOC_DH(2*m*sizeof(MPI_Request)); CHECK_V_ERROR;
    rcv_status = (MPI_Status*)MALLOC_DH(2*m*sizeof(MPI_Status)); CHECK_V_ERROR;

    for (i=0; i<m; ++i) {

      /* error check for empty row */
      int count = rp[i+1] - rp[i];
      if (! count) {
        sprintf(msgBuf_dh, "local row %i of %i is empty!", i+1, m);
        SET_V_ERROR(msgBuf_dh);
      }

      MPI_Irecv(cval+rp[i], count, MPI_INT, 0, CVAL_TAG, comm_dh, rcv_req+2*i);
      MPI_Irecv(aval+rp[i], count, MPI_DOUBLE, 0, AVAL_TAG, comm_dh, rcv_req+2*i+1);
    }
  }

  /* wait for all sends/receives to finish */
  if (myid_dh == 0) {
    MPI_Waitall(m*2, send_req, send_status);
  }
  MPI_Waitall(2*B->m, rcv_req, rcv_status);

  /* clean up */
  if (rowLengths != NULL) { FREE_DH(rowLengths); CHECK_V_ERROR; }
  if (o2n_row != NULL) { FREE_DH(o2n_row); CHECK_V_ERROR; }
  if (n2o_col != NULL) { FREE_DH(n2o_col); CHECK_V_ERROR; }
  if (rowToBlock != NULL) {FREE_DH(rowToBlock); CHECK_V_ERROR; }
  if (send_req != NULL) { FREE_DH(send_req); CHECK_V_ERROR; }
  if (rcv_req != NULL) { FREE_DH(rcv_req); CHECK_V_ERROR; }
  if (send_status != NULL) { FREE_DH(send_status); CHECK_V_ERROR; }
  if (rcv_status != NULL) { FREE_DH(rcv_status); CHECK_V_ERROR; }
  if (beg_row != NULL) { FREE_DH(beg_row); CHECK_V_ERROR; }
  if (row_count != NULL) { FREE_DH(row_count); CHECK_V_ERROR; }
  if (C != NULL) { Mat_dhDestroy(C); ERRCHKA; }

  *Bout = B;

  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "partition_and_distribute_private"
void partition_and_distribute_private(Mat_dh A, Mat_dh *Bout)
{
  START_FUNC_DH
  Mat_dh B = NULL;
  int i, m;
  int *rowLengths = NULL;
  int *o2n_row = NULL, *n2o_col = NULL, *rowToBlock = NULL;
  MPI_Request *send_req = NULL;
  MPI_Request *rcv_req = NULL;
  MPI_Status  *send_status = NULL;
  MPI_Status  *rcv_status = NULL;

  MPI_Barrier(comm_dh);

  /* broadcast number of rows to all processors */
  if (myid_dh == 0)  m = A->m;
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* broadcast number of nonzeros in each row to all processors */
  rowLengths = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  if (myid_dh == 0) {
    int *tmp = A->rp; 
    for (i=0; i<m; ++i) {
      rowLengths[i] = tmp[i+1] - tmp[i];
    }
  }
  MPI_Bcast(rowLengths, m, MPI_INT, 0, comm_dh);

  /* partition matrix */
  rowToBlock = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;

  if (myid_dh == 0) {
    o2n_row = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
    mat_partition_private(A, np_dh, o2n_row, rowToBlock); CHECK_V_ERROR;
  }

  /* broadcast partitiioning information to all processors */
  MPI_Bcast(rowToBlock, m, MPI_INT, 0, comm_dh);

  /* allocate storage for local portion of matrix */
  mat_par_read_allocate_private(&B, m, rowLengths, rowToBlock); CHECK_V_ERROR;

  /* root sends each processor its portion of the matrix */
  if (myid_dh == 0) {
    int *cval = A->cval, *rp = A->rp;
    double *aval = A->aval;
    send_req = (MPI_Request*)MALLOC_DH(2*m*sizeof(MPI_Request)); CHECK_V_ERROR;
    send_status = (MPI_Status*)MALLOC_DH(2*m*sizeof(MPI_Status)); CHECK_V_ERROR;
    for (i=0; i<m; ++i) {
      int owner = rowToBlock[i];
      int count = rp[i+1]-rp[i];

      /* error check for empty row */
      if (! count) {
        sprintf(msgBuf_dh, "row %i of %i is empty!", i+1, m);
        SET_V_ERROR(msgBuf_dh);
      }

      MPI_Isend(cval+rp[i], count, MPI_INT, owner, CVAL_TAG, comm_dh, send_req+2*i);
      MPI_Isend(aval+rp[i], count, MPI_DOUBLE, owner, AVAL_TAG, comm_dh, send_req+2*i+1);
    }
  } 

  /* all processors receive their local rows */
  { int *cval = B->cval;
    int *rp = B->rp;
    double *aval = B->aval;
    m = B->m;

    rcv_req = (MPI_Request*)MALLOC_DH(2*m*sizeof(MPI_Request)); CHECK_V_ERROR;
    rcv_status = (MPI_Status*)MALLOC_DH(2*m*sizeof(MPI_Status)); CHECK_V_ERROR;

    for (i=0; i<m; ++i) {

      /* error check for empty row */
      int count = rp[i+1] - rp[i];
      if (! count) {
        sprintf(msgBuf_dh, "local row %i of %i is empty!", i+1, m);
        SET_V_ERROR(msgBuf_dh);
      }

      MPI_Irecv(cval+rp[i], count, MPI_INT, 0, CVAL_TAG, comm_dh, rcv_req+2*i);
      MPI_Irecv(aval+rp[i], count, MPI_DOUBLE, 0, AVAL_TAG, comm_dh, rcv_req+2*i+1);
    }
  }

  /* wait for all sends/receives to finish */
  if (myid_dh == 0) {
    MPI_Waitall(m*2, send_req, send_status);
  }
  MPI_Waitall(2*B->m, rcv_req, rcv_status);

  /* clean up */
  if (rowLengths != NULL) { FREE_DH(rowLengths); CHECK_V_ERROR; }
  if (o2n_row != NULL) { FREE_DH(o2n_row); CHECK_V_ERROR; }
  if (n2o_col != NULL) { FREE_DH(n2o_col); CHECK_V_ERROR; }
  if (rowToBlock != NULL) {FREE_DH(rowToBlock); CHECK_V_ERROR; }
  if (send_req != NULL) { FREE_DH(send_req); CHECK_V_ERROR; }
  if (rcv_req != NULL) { FREE_DH(rcv_req); CHECK_V_ERROR; }
  if (send_status != NULL) { FREE_DH(send_status); CHECK_V_ERROR; }
  if (rcv_status != NULL) { FREE_DH(rcv_status); CHECK_V_ERROR; }

  *Bout = B;

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "mat_par_read_allocate_private"
void mat_par_read_allocate_private(Mat_dh *Aout, int n, int *rowLengths, int *rowToBlock)
{
  START_FUNC_DH
  Mat_dh A;
  int i, m, nz, beg_row, *rp, idx;

  Mat_dhCreate(&A); CHECK_V_ERROR;
  *Aout =  A;
  A->n = n;

  /* count number of rows owned by this processor */
  m = 0;
  for (i=0; i<n; ++i) {
    if (rowToBlock[i] == myid_dh) ++m;
  }
  A->m = m;

  /* compute global numbering of first  locally owned row */
  beg_row = 0;
  for (i=0; i<n; ++i) {
    if (rowToBlock[i] < myid_dh) ++beg_row;
  }
  A->beg_row = beg_row;

  /* allocate storage for row-pointer array */
  A->rp = rp = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  rp[0] = 0;

  /* count number of nonzeros owned by this processor, and form rp array */
  nz = 0;
  idx = 1;
  for (i=0; i<n; ++i) {
    if (rowToBlock[i] == myid_dh) {
      nz += rowLengths[i];
      rp[idx++] = nz;
    }
  }

  /* allocate storage for column indices and values arrays */
  A->cval = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  A->aval = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "mat_partition_private"
void mat_partition_private(Mat_dh A, int blocks, int *o2n_row, int *rowToBlock)
{
  START_FUNC_DH
  int i, j, n = A->n;
  int rpb = n/blocks;   /* rows per block (except possibly last block) */
  int idx = 0;

  while (rpb*blocks < n) ++rpb;

  if (rpb*(blocks-1) == n) {
    --rpb;
    printf_dh("adjusted rpb to: %i\n", rpb);
  }

  for (i=0; i<n; ++i) o2n_row[i] = i;

  /* assign all rows to blocks, except for last block, which may
     contain less than "rpb" rows
   */
  for (i=0; i<blocks-1; ++i) {
    for (j=0; j<rpb; ++j) {
      rowToBlock[idx++] = i;
    }  
  }
 
  /* now deal with the last block in the partition */
  i = blocks - 1;
  while (idx < n) rowToBlock[idx++] = i;

  END_FUNC_DH
}


/* may produce incorrect result if input is not triangular! */
#undef __FUNC__
#define __FUNC__ "make_full_private"
void make_full_private(int m, int **rpIN, int **cvalIN, double **avalIN)
{
  START_FUNC_DH
  int i, j, *rpNew, *cvalNew, *rp = *rpIN, *cval = *cvalIN;
  double *avalNew, *aval = *avalIN;
  int nz, *rowCounts = NULL;

  /* count the number of nonzeros in each row */
  rowCounts = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<=m; ++i) rowCounts[i] = 0;

  for (i=0; i<m; ++i) {
    for (j=rp[i]; j<rp[i+1]; ++j) {
      int col = cval[j];
      rowCounts[i+1] += 1;
      if (col != i) rowCounts[col+1] += 1;
    }
  }

  /* prefix sum to form row pointers for full representation */
  rpNew = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  for (i=1; i<=m; ++i) rowCounts[i] += rowCounts[i-1];
  memcpy(rpNew, rowCounts, (m+1)*sizeof(int));

  /* form full representation */
  nz = rpNew[m];

  cvalNew = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  avalNew = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) {
    for (j=rp[i]; j<rp[i+1]; ++j) {
      int col = cval[j];
      double val  = aval[j];

      cvalNew[rowCounts[i]] = col;
      avalNew[rowCounts[i]] = val;
      rowCounts[i] += 1;
      if (col != i) {
        cvalNew[rowCounts[col]] = i;
        avalNew[rowCounts[col]] = val;
        rowCounts[col] += 1;
      }
    }
  }

  if (rowCounts != NULL) { FREE_DH(rowCounts); CHECK_V_ERROR; }
  FREE_DH(cval); CHECK_V_ERROR;
  FREE_DH(rp); CHECK_V_ERROR;
  FREE_DH(aval); CHECK_V_ERROR;
  *rpIN = rpNew;
  *cvalIN = cvalNew;
  *avalIN = avalNew;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "make_symmetric_private"
void make_symmetric_private(int m, int **rpIN, int **cvalIN, double **avalIN)
{
  START_FUNC_DH
  int i, j, *rpNew, *cvalNew, *rp = *rpIN, *cval = *cvalIN;
  double *avalNew, *aval = *avalIN;
  int nz, *rowCounts = NULL;
  int *rpTrans, *cvalTrans;
  int *work;
  double *avalTrans;
  int nzCount = 0, transCount = 0;

  mat_dh_transpose_private(m, rp, &rpTrans,
                           cval, &cvalTrans, aval, &avalTrans); CHECK_V_ERROR;

  /* count the number of nonzeros in each row */
  work = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) work[i] = -1;
  rowCounts = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<=m; ++i) rowCounts[i] = 0;

  for (i=0; i<m; ++i) {
    int ct = 0;
    for (j=rp[i]; j<rp[i+1]; ++j) {
      int col = cval[j];
      work[col] = i;
      ++ct;
      ++nzCount;
    }
    for (j=rpTrans[i]; j<rpTrans[i+1]; ++j) {
      int col = cvalTrans[j];
      if (work[col] != i) {
        ++ct;
        ++transCount;
      }
    }
    rowCounts[i+1] = ct;
  }

  /*---------------------------------------------------------
   * if matrix is already symmetric, do nothing
   *---------------------------------------------------------*/
  if (transCount == 0) {
    printf("make_symmetric_private: matrix is already structurally symmetric!\n");
    FREE_DH(rpTrans); CHECK_V_ERROR;
    FREE_DH(cvalTrans); CHECK_V_ERROR;
    FREE_DH(avalTrans); CHECK_V_ERROR;
    FREE_DH(work); CHECK_V_ERROR;
    FREE_DH(rowCounts); CHECK_V_ERROR;
    goto END_OF_FUNCTION;
  } 

  /*---------------------------------------------------------
   * otherwise, finish symmetrizing
   *---------------------------------------------------------*/
    else {
    printf("original nz= %i\n", rp[m]);
    printf("zeros added= %i\n", transCount);
    printf("ratio of added zeros to nonzeros = %0.2f (assumes all original entries were nonzero!)\n", 
                 (double)transCount/(double)(nzCount) );
  }

  /* prefix sum to form row pointers for full representation */
  rpNew = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  for (i=1; i<=m; ++i) rowCounts[i] += rowCounts[i-1];
  memcpy(rpNew, rowCounts, (m+1)*sizeof(int));
  for (i=0; i<m; ++i) work[i] = -1;

  /* form full representation */
  nz = rpNew[m];
  cvalNew = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  avalNew = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) work[i] = -1;

  for (i=0; i<m; ++i) {
    for (j=rp[i]; j<rp[i+1]; ++j) {
      int col = cval[j];
      double val  = aval[j];
      work[col] = i;
      cvalNew[rowCounts[i]] = col;
      avalNew[rowCounts[i]] = val;
      rowCounts[i] += 1;
    }
    for (j=rpTrans[i]; j<rpTrans[i+1]; ++j) {
      int col = cvalTrans[j];
      if (work[col] != i) {
        cvalNew[rowCounts[i]] = col;
        avalNew[rowCounts[i]] = 0.0;
        rowCounts[i] += 1;
      }
    }
  }

  if (rowCounts != NULL) { FREE_DH(rowCounts); CHECK_V_ERROR; }
  FREE_DH(work); CHECK_V_ERROR;
  FREE_DH(cval); CHECK_V_ERROR;
  FREE_DH(rp); CHECK_V_ERROR;
  FREE_DH(aval); CHECK_V_ERROR;
  FREE_DH(cvalTrans); CHECK_V_ERROR;
  FREE_DH(rpTrans); CHECK_V_ERROR;
  FREE_DH(avalTrans); CHECK_V_ERROR;
  *rpIN = rpNew;
  *cvalIN = cvalNew;
  *avalIN = avalNew;

END_OF_FUNCTION: ;

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "profileMat"
void profileMat(Mat_dh A)
{
  START_FUNC_DH
  Mat_dh B = NULL;
  int type;
  int m;
  int i, j;
  int *work1;
  double *work2;
  bool isStructurallySymmetric = true;
  bool isNumericallySymmetric = true;
  bool is_Triangular = false;
  int zeroCount = 0, nz;

  if (myid_dh > 0) {
    SET_V_ERROR("only for a single MPI task!");
  }

  m = A->m;

  printf("\nYY----------------------------------------------------\n");

  /* count number of explicit zeros */
  nz = A->rp[m];
  for (i=0; i<nz; ++i) {
    if (A->aval[i] == 0) ++zeroCount;
  }
  printf("YY  row count:      %i\n", m);
  printf("YY  nz count:       %i\n", nz);
  printf("YY  explicit zeros: %i (entire matrix)\n", zeroCount);

  /* count number of missing or zero diagonals */
  { int m_diag = 0, z_diag = 0;
    for (i=0; i<m; ++i) {
      bool flag = true;
      for (j=A->rp[i]; j<A->rp[i+1]; ++j) {
        int col = A->cval[j];

        /* row has an explicit diagonal element */
        if (col == i) {          
          double val = A->aval[j];
          flag = false;
          if (val == 0.0) ++z_diag;
          break;
        }
      }

      /* row has an implicit zero diagonal element */
      if (flag) ++m_diag;
    }
    printf("YY  missing diagonals:   %i\n", m_diag);
    printf("YY  explicit zero diags: %i\n", z_diag);
  }

  /* check to see if matrix is triangular */
  type = isTriangular(m, A->rp, A->cval); CHECK_V_ERROR;
  if (type == IS_UPPER_TRI) {
    printf("YY  matrix is upper triangular\n");
    is_Triangular = true;
    goto END_OF_FUNCTION;
  } else if (type == IS_LOWER_TRI) {
    printf("YY  matrix is lower triangular\n");
    is_Triangular = true;
    goto END_OF_FUNCTION;
  }

  /* if not triangular, count nz in each triangle */
  { int unz = 0, lnz = 0;
    for (i=0; i<m; ++i) {
      for (j=A->rp[i]; j<A->rp[i+1]; ++j) {
        int col = A->cval[j];
        if (col < i) ++lnz;
        if (col > i) ++unz;
      }
    }
    printf("YY  strict upper triangular nonzeros: %i\n", unz);
    printf("YY  strict lower triangular nonzeros: %i\n", lnz);
  }
 
   
  

  Mat_dhTranspose(A, &B); CHECK_V_ERROR;

  /* check for structural and numerical symmetry */

  work1 = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  work2 = (double*)MALLOC_DH(m*sizeof(double)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) work1[i] = -1;
  for (i=0; i<m; ++i) work2[i] = 0.0;

  for (i=0; i<m; ++i) {
    for (j=A->rp[i]; j<A->rp[i+1]; ++j) {
      int col = A->cval[j];
      double val = A->aval[j];
      work1[col] = i;
      work2[col] = val;
    }
    for (j=B->rp[i]; j<B->rp[i+1]; ++j) {
      int col = B->cval[j];
      double val = B->aval[j];

      if (work1[col] != i) {
        isStructurallySymmetric = false;
        isNumericallySymmetric = false;
        goto END_OF_FUNCTION;
      }
      if (work2[col] != val) {
        isNumericallySymmetric = false;
        work2[col] = 0.0;
      }
    }
  }


END_OF_FUNCTION: ;

  if (! is_Triangular) {
    printf("YY  matrix is NOT triangular\n");
    if (isStructurallySymmetric) {
      printf("YY  matrix IS structurally symmetric\n");
    } else {
      printf("YY  matrix is NOT structurally symmetric\n");
    }
    if (isNumericallySymmetric) {
      printf("YY  matrix IS numerically symmetric\n");
    } else {
      printf("YY  matrix is NOT numerically symmetric\n");
    }
  }

  if (work1 != NULL) { FREE_DH(work1); CHECK_V_ERROR; }
  if (work2 != NULL) { FREE_DH(work2); CHECK_V_ERROR; }
  if (B != NULL) { Mat_dhDestroy(B); CHECK_V_ERROR; }

  printf("YY----------------------------------------------------\n");

  END_FUNC_DH
}
