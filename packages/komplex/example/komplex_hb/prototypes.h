/*@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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

void read_hb(char *data_file, int *proc_config,
	      int *N_global, int *n_nonzeros, 
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);

void read_coo(char *data_file, int *proc_config,
	      int *N_global, int *n_nonzeros,
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);


void distrib_msr_matrix(int *proc_config,
	      int N_global, int *n_nonzeros, 
           int *N_update, int **update, 
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);

void distrib_vbr_matrix(int *proc_config,
	      int N_global, int N_blk_global, 
           int *n_nonzeros,  int *n_blk_nonzeros,
           int *N_update, int **update, 
	      double **val, int **indx, int **rpntr, int **cpntr,
           int **bpntr, int **bindx,
	      double **x, double **b, double **xexact);

void create_vbr(char *part_file, int *proc_config, 
                int *N_global, int *N_blk_global, 
                int *n_nonzeros, int *n_blk_nonzeros,
                int *N_update, int **update,
		      int *bindx_msr, double *val_msr,
		      double **val, int **indx, int **rpntr, int **cpntr,
		      int **bpntr, int **bindx);

double smsrres (int m, int n, 
	      double *val, int *indx, 
	      double *xlocal, double *x, double *b);

double scscres (int isym, int m, int n, 
	      double *val, int *indx, int *pntr,
	      double *x, double *b);

void  scscmv (int isym, int m, int n, 
	      double *val, int *indx, int *pntr,
	      double *x, double *b);

double svbrres (int m, int n, int m_blk,
		double *val, int *indx, int *bindx, int *rpntr,
		int *cpntr, int *bpntrb, int *bpntre,
		double *x, double *b);
