/***************************************************************************
                                decsol.c
                             -------------------
    modified for C by:   : Blake Ashby
    last modified        : Nov 15, 2002
    email                : bmashby@stanford.edu

Helper routines for StiffIntegratorT class. Modified from original Fortran
code provided by:

	E. Hairer & G. Wanner
	Universite de Geneve, dept. de Mathematiques
	CH-1211 GENEVE 4, SWITZERLAND
	E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

 ***************************************************************************/

#include "decsol.h"
#include "math.h"

int dec(const int n, double **A, int *ip)
{

/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination

	Input:
		n		order of matrix
		A		matrix to be triangularized

	Output:
		A[i][j], i <= j		upper triangular factor, U
		A[i][j], i > j		multipliers = lower triangular factor, i - l
		ip[k], k < n - 1	index of k-th pivot rowf
		ip[n-1]				(-1)^(number of interchanges) or 0
		ier 				0 if matrix A is nonsingular, or k if found
								to be singular at stage k

	Use sol to obtain solution of linear system

	determ(A) = ip[n-1]*A[0][0]*A[1][1]*...*A[n-1][n-1]

	If ip[n-1] = 0, A is singular, sol will divide by zero

	Reference:
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/

	int kp1, m, nm1, k, ier;
	double t;

	ier = 0;
	ip[n-1] = 1;
	if (n != 1) {
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = k;
			for (int i = kp1; i < n; i++) {
				if (fabs(A[i][k]) > fabs(A[m][k])) m = i;
			}
			ip[k] = m;
			t = A[m][k];
			if (m != k) {
				ip[n-1] = -ip[n-1];
				A[m][k] = A[k][k];
				A[k][k] = t;
			}
			if (t == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			t = 1.0/t;
			for (int i = kp1; i < n; i++) A[i][k] *= -t;
			for (int j = kp1; j < n; j++) {
				t = A[m][j];
				A[m][j] = A[k][j];
				A[k][j] = t;
				if (t != 0.0)
					for (int i = kp1; i < n; i++)
						A[i][j] += A[i][k]*t;
			}
		}
	}
	k = n;
	if (A[n-1][n-1] == 0.0) {
		ier = k;
		ip[n-1] = 0;
	}

	return (ier);

} //dec

void sol(const int n, double **A, double *b, int *ip)
{

/*-----------------------------------------------------------------------
	Solution of linear system A*x = b

	Input:
		n		order of matrix
		A		triangularized matrix obtained from dec
		b		right hand side vector
		ip		pivot vector obtained from dec

	Do not use if dec has set ier != 0

	Output:
		b		solution vector, x

-----------------------------------------------------------------------*/

	int m, kb, km1, nm1, kp1;
	double t;

	if (n != 1) {
		nm1 = n - 1;
		for (int k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = ip[k];
			t = b[m];
			b[m] = b[k];
			b[k] = t;
			for (int i = kp1; i < n; i++)
				b[i] += A[i][k]*t;
		}
		for (int k = 0; k < nm1; k++) {
			km1 = n - k - 2;
			kb = km1 + 1;
			b[kb] = b[kb]/A[kb][kb];
			t = -b[kb];
			for (int i = 0; i <= km1; i++) 
				b[i] += A[i][kb]*t;
		}
	}
	b[0] = b[0]/A[0][0];
	
	return;
	
} //sol


int dech(const int n, double **A, int lb, int *ip)
{

/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination of a Hessenberg
		matrix with lower bandwidth lb
	
	Input:
		n		order of matrix
		A		matrix to be triangularized
		lb		lower bandwidth of A (diagonal is not counted, lb >= 1)
		
	Output:
		A[i][j], i <= j		upper triangular factor, U
		A[i][j], i > j		multipliers = lower triangular factor, i - l
		ip[k], k < n - 1	index of k-th pivot row
		ip[n-1] 			(-1)^(number of interchanges) or 0
		ier 				0 if matrix A is nonsingular, or k if found 
								to be singular at stage k

	Use solh to obtain solution of linear system
	
	determ(A) = ip[n-1]*A[0][0]*A[1][1]*...*A[n-1][n-1]
	
	If ip[n-1] = 0, A is singular, solh will divide by zero
	
	Reference:
		This is a slight modification of
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/	  

	int kp1, m, nm1, k, na, ier;
	double t;

	ier = 0;
	ip[n-1] = 1;
	if (n != 1) {
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = k;
			na = min(n, lb+k+1);
			for (int i = kp1; i < na; i++) {
				if (fabs(A[i][k]) > fabs(A[m][k])) m = i;
			}
			ip[k] = m;
			t = A[m][k];
			if (m != k) {
				ip[n-1] = -ip[n-1];
				A[m][k] = A[k][k];
				A[k][k] = t;
			}
			if (t == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			t = 1.0/t;
			for (int i = kp1; i < n; i++) A[i][k] *= -t;
			for (int j = kp1; j < n; j++) {
				t = A[m][j];
				A[m][j] = A[k][j];
				A[k][j] = t;
				if (t != 0.0)
					for (int i = kp1; i < na; i++) 
						A[i][j] += A[i][k]*t;	
			}	
		}
	}
	k = n;
	if (A[n-1][n-1] == 0.0) {
		ier = k;
		ip[n-1] = 0;
	}
	
	return (ier);

} //dech


void solh(const int n, double **A, int lb, double *b, int *ip)
{

/*-----------------------------------------------------------------------
	Solution of linear system A*x = b -- Hessenberg matrix
	
	Input:
		n		order of matrix
		A		triangularized matrix obtained from dech
		lb		lower bandwidth of A
		b		right hand side vector
		ip		pivot vector obtained from dec
		
	Do not use if dec has set ier != 0
	
	Output:
		b		solution vector, x

-----------------------------------------------------------------------*/

	int m, kb, km1, nm1, kp1, na;
	double t;

	if (n != 1) {
		nm1 = n - 1;
		for (int k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = ip[k];
			t = b[m];
			b[m] = b[k];
			b[k] = t;
			na = min(n, lb+k+1);
			for (int i = kp1; i < na; i++)
				b[i] += A[i][k]*t;
		}
		for (int k = 0; k < nm1; k++) {
			km1 = n - k - 2;
			kb = km1 + 1;
			b[kb] = b[kb]/A[kb][kb];
			t = -b[kb];
			for (int i = 0; i <= km1; i++) 
				b[i] += A[i][kb]*t;
		}
	}
	b[0] = b[0]/A[0][0];
	
	return;

} //solh

int decc(const int n, double **AR, double **AI, int *ip)
{

/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination
	------ Modification for complex matrices -------
	
	Input:
		n			order of matrix
		AR, AI		matrix to be triangularized
		
	Output:
		AR[i][j], i <= j	upper triangular factor, U; real part
		AI[i][j], i <= j	upper triangular factor, U; imaginary part
		AR[i][j], i > j		multipliers = lower triangular factor, i - l 
								real part
		AI[i][j], i > j		multipliers = lower triangular factor, i - l
								imaginary part
		ip[k], k < n - 1	index of k-th pivot row
		ip[n-1]				(-1)^(number of interchanges) or 0
		ier					0 if matrix A is nonsingular, or k if found 
								to be singular at stage k
			
	Use solc to obtain solution of linear system
	
	If ip[n-1] = 0, A is singular, sol will divide by zero
	
	Reference:
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/	  

	int kp1, m, nm1, k, ier;
	double tr, ti, den, prodr, prodi;

	ier = 0;
	ip[n-1] = 1;
	if (n != 1) {
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = k;
			for (int i = kp1; i < n; i++) {
				if ((fabs(AR[i][k]) + fabs(AI[i][k])) >
					(fabs(AR[m][k]) + fabs(AI[m][k]))) 
					m = i;
			}
			ip[k] = m;
			tr = AR[m][k];
			ti = AI[m][k];
			if (m != k) {
				ip[n-1] = -ip[n-1];
				AR[m][k] = AR[k][k];
				AI[m][k] = AI[k][k];
				AR[k][k] = tr;
				AI[k][k] = ti;
			}
			if ((fabs(tr) + fabs(ti)) == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			den = tr*tr + ti*ti;
			tr = tr/den;
			ti = -ti/den;
			for (int i = kp1; i < n; i++) {
				prodr = AR[i][k]*tr - AI[i][k]*ti;
				prodi = AI[i][k]*tr + AR[i][k]*ti;
				AR[i][k] = -prodr;
				AI[i][k] = -prodi;
			}
			for (int j = kp1; j < n; j++) {
				tr = AR[m][j];
				ti = AI[m][j];
				AR[m][j] = AR[k][j];
				AI[m][j] = AI[k][j];
				AR[k][j] = tr;
				AI[k][j] = ti;
				if ((fabs(tr) + fabs(ti)) == 0.0) {
				}
				else if (ti == 0.0) {
					for (int i = kp1; i < n; i++) {
						prodr = AR[i][k]*tr;
						prodi = AI[i][k]*tr;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
				else if (tr == 0.0) {
					for (int i = kp1; i < n; i++) {
						prodr = -AI[i][k]*ti;
						prodi = AR[i][k]*ti;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
				else {
					for (int i = kp1; i < n; i++) {
						prodr = AR[i][k]*tr - AI[i][k]*ti;
						prodi = AI[i][k]*tr + AR[i][k]*ti;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
			}
		}	
	}	
	k = n;
	if ((fabs(AR[n-1][n-1]) + fabs(AI[n-1][n-1])) == 0.0) {
		ier = k;
		ip[n-1] = 0;
	}

	return (ier);

} //decc

void solc(const int n, double **AR, double **AI, double *br,
	double *bi, int *ip)
{

/*-----------------------------------------------------------------------
	Solution of linear system A*x = b
	
	Input:
		n			order of matrix
		AR, AI		triangularized matrix obtained from decc
		br, bi		right hand side vector
		ip			pivot vector obtained from dec
		
	Do not use if decc has set ier != 0
	
	Output:
		br, bi		solution vector, x

-----------------------------------------------------------------------*/

	int nm1, kp1, m, kb, km1;
	double den, prodr, prodi, tr, ti;
	  
	if (n != 1) {
		nm1 = n - 1;
		for (int k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = ip[k];
			tr = br[m];
			ti = bi[m];
			br[m] = br[k];
			bi[m] = bi[k];
			br[k] = tr;
			bi[k] = ti;
			for (int i = kp1; i < n; i++) {
				prodr = AR[i][k]*tr - AI[i][k]*ti;
				prodi = AI[i][k]*tr + AR[i][k]*ti;
				br[i] += prodr;
				bi[i] += prodi;
			}
		} 
 		for (int k = 0; k < nm1; k++) {
			km1 = n - k - 2;
			kb = km1 + 1;
			den = AR[kb][kb]*AR[kb][kb] + AI[kb][kb]*AI[kb][kb];
			prodr = br[kb]*AR[kb][kb] + bi[kb]*AI[kb][kb];
			prodi = bi[kb]*AR[kb][kb] - br[kb]*AI[kb][kb];
			br[kb] = prodr/den;
			bi[kb] = prodi/den;
			tr = -br[kb];
			ti = -bi[kb];
			for (int i = 0; i <= km1; i++) {
				prodr = AR[i][kb]*tr - AI[i][kb]*ti;
				prodi = AI[i][kb]*tr + AR[i][kb]*ti;
				br[i] += prodr;
				bi[i] += prodi;
			}
		}
	}
	den = AR[0][0]*AR[0][0] + AI[0][0]*AI[0][0];
	prodr = br[0]*AR[0][0] + bi[0]*AI[0][0];
	prodi = bi[0]*AR[0][0] - br[0]*AI[0][0];
	br[0] = prodr/den;
	bi[0] = prodi/den;
	
	return;
		
} // solc

int dechc(const int n, double **AR, double **AI, int lb, int *ip)
{

/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination
	------ Modification for complex matrices -------
	
	Input:
		n			order of matrix
		AR, AI		matrix to be triangularized
		
	Output:
		AR[i][j], i <= j	upper triangular factor, U; real part
		AI[i][j], i <= j	upper triangular factor, U; imaginary part
		AR[i][j], i > j		multipliers = lower triangular factor, i - l 
								real part
		AI[i][j], i > j		multipliers = lower triangular factor, i - l 
								imaginary part
		lb					lower bandwidth of A (diagonal not counted), lb >= 1
		ip[k], k < n - 1	index of k-th pivot row
		ip[n-1]				(-1)^(number of interchanges) or 0
		ier					0 if matrix A is nonsingular, or k if found
								to be singular at stage k
			
	Use solhc to obtain solution of linear system

	If ip[n-1] = 0, A is singular, solhc will divide by zero

	Reference:
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/	  

	int kp1, m, nm1, k, na, ier;
	double tr, ti, den, prodr, prodi;
 
	ier = 0;
	ip[n-1] = 1;
	if ((n != 1) && (lb != 0)) {
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = k;
			na = min(n, lb+k+1);
			for (int i = kp1; i < na; i++) {
				if ((fabs(AR[i][k]) + fabs(AI[i][k])) > 
					(fabs(AR[m][k]) + fabs(AI[m][k]))) 
					m = i;
			}
			ip[k] = m;
			tr = AR[m][k];
			ti = AI[m][k];
			if (m != k) {
				ip[n-1] = -ip[n-1];
				AR[m][k] = AR[k][k];
				AI[m][k] = AI[k][k];
				AR[k][k] = tr;
				AI[k][k] = ti;
			}
			if ((fabs(tr) + fabs(ti)) == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			den = tr*tr + ti*ti;
			tr = tr/den;
			ti = -ti/den;
			for (int i = kp1; i < na; i++) {
				prodr = AR[i][k]*tr - AI[i][k]*ti;
				prodi = AI[i][k]*tr + AR[i][k]*ti;
				AR[i][k] = -prodr;
				AI[i][k] = -prodi;
			}
			for (int j = kp1; j < n; j++) {
				tr = AR[m][j];
				ti = AI[m][j];
				AR[m][j] = AR[k][j];
				AI[m][j] = AI[k][j];
				AR[k][j] = tr;
				AI[k][j] = ti;
				if ((fabs(tr) + fabs(ti)) == 0.0) {
				}
				else if (ti == 0.0) {
					for (int i = kp1; i < na; i++) {
						prodr = AR[i][k]*tr;
						prodi = AI[i][k]*tr;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
				else if (tr == 0.0) {
					for (int i = kp1; i < na; i++) {
						prodr = -AI[i][k]*ti;
						prodi = AR[i][k]*ti;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
				else {
					for (int i = kp1; i < na; i++) {
						prodr = AR[i][k]*tr - AI[i][k]*ti;
						prodi = AI[i][k]*tr + AR[i][k]*ti;
						AR[i][j] += prodr;
						AI[i][j] += prodi;
					}
				}
			}
		}	
	}	
	k = n;
	if ((fabs(AR[n-1][n-1]) + fabs(AI[n-1][n-1])) == 0.0) {
		ier = k;
		ip[n-1] = 0;
	}

	return (ier);

} // dechc

void solhc(const int n, double **AR, double **AI, int lb,
	double *br, double *bi, int *ip)
{

/*-----------------------------------------------------------------------
	Solution of linear system A*x = b

	Input:
		n			order of matrix
		AR, AI		triangularized matrix obtained from dec
		br, bi		right hand side vector
		lb			lower bandwidth of A
		ip			pivot vector obtained from dec
		
	Do not use if dechc has set ier != 0
	
	Output:
		br, bi		solution vector, x

-----------------------------------------------------------------------*/

	int nm1, kp1, m, kb, km1;
	double den, prodr, prodi, tr, ti;
	  
	if (n != 1) {
		nm1 = n - 1;
		if (lb != 0) {
			for (int k = 0; k < nm1; k++) {
				kp1 = k + 1;
				m = ip[k];
				tr = br[m];
				ti = bi[m];
				br[m] = br[k];
				bi[m] = bi[k];
				br[k] = tr;
				bi[k] = ti;
				for (int i = kp1; i < min(n, lb+k+1); i++) {
					prodr = AR[i][k]*tr - AI[i][k]*ti;
					prodi = AI[i][k]*tr + AR[i][k]*ti;
					br[i] += prodr;
					bi[i] += prodi;
				}
			}
		}
 		for (int k = 0; k < nm1; k++) {
			km1 = n - k - 2;
			kb = km1 + 1;
			den = AR[kb][kb]*AR[kb][kb] + AI[kb][kb]*AI[kb][kb];
			prodr = br[kb]*AR[kb][kb] + bi[kb]*AI[kb][kb];
			prodi = bi[kb]*AR[kb][kb] - br[kb]*AI[kb][kb];
			br[kb] = prodr/den;
			bi[kb] = prodi/den;
			tr = -br[kb];
			ti = -bi[kb];
			for (int i = 0; i < km1; i++) {
				prodr = AR[i][kb]*tr - AI[i][kb]*ti;
				prodi = AI[i][kb]*tr + AR[i][kb]*ti;
				br[i] += prodr;
				bi[i] += prodi;
			}
		}
	}
	den = AR[0][0]*AR[0][0] + AI[0][0]*AI[0][0];
	prodr = br[0]*AR[0][0] + bi[0]*AI[0][0];
	prodi = bi[0]*AR[0][0] - br[0]*AI[0][0];
	br[0] = prodr/den;
	bi[0] = prodi/den;
	
	return;

}

int decb(const int n, double **A, int ml, int mu, int *ip)
{

/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination of a banded
		matrix with lower bandwidth ml and upper bandwidth mu
	
	Input:
		n		order of matrix
		A		contains the matrix in band storage.
				The columns of the matrix are stored in the columns
				of A and the diagonals of the matrix are stored in 
				rows ml	through 2*ml + mu of A.
		ml		lower bandwidth of A (diagonal is not counted)
		mu		upper bandwidth of A (diagonal is not counted)
				
	Output:
		A 		upper triangular matrix in band storage and the
				multipliers which were used to obtain it
		ip		index vector of pivot indices
		ip[n-1]	(-1)^(number of interchanges) or 0
		ier		0 if matrix A is nonsingular, or k if found to be 
					singular at stage k
			
	Use solb to obtain solution of linear system
	
	determ(A) = ip[n-1]*A[md][0]*A[md][1]*...*A[md][n-1] with
			md = ml + mu
	
	If ip[n-1] = 0, A is singular, solb will divide by zero
	
	Reference:
		This is a modification of:
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/

	int md, md1, mdl, ju, nm1, k, kp1, m, mm, jk, ijk, ier;
	double t;

	ier = 0;
	ip[n-1] = 1;
	md = ml + mu;
	md1 = md + 1;
	ju = 0;
	if ((n != 1) && (ml != 0)) {
		if (n >= mu+2)
			for (int j = mu + 1; j < n; j++) 
				for (int i = 0; i < ml; i++) 
					A[i][j] = 0.0;
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = md;
			mdl = min(ml, n - k - 1) + md;   
			for (int i = md1; i <= mdl; i++) { 
				if (fabs(A[i][k]) > fabs(A[m][k])) m = i;
			}
			ip[k] = m + k - md;
			t = A[m][k];
			if (m != md) {
				ip[n-1] = -ip[n-1];
				A[m][k] = A[md][k];
				A[md][k] = t;
			}
			if (t == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			t = 1.0/t;
			for (int i = md1; i <= mdl; i++) A[i][k] *= -t; 
			ju = min(max(ju, mu+ip[k]+1), n);
			mm = md;
			if (ju >= kp1) {
				for (int j = kp1; j < ju; j++) {
					m = m - 1;
					mm = mm - 1;
					t = A[m][j];
					if (m != mm) {
						A[m][j] = A[mm][j];
						A[mm][j] = t;
					}
					if (t != 0.0) {
						jk = j - k;
						for (int i = md1; i <= mdl; i++) {
							ijk = i - jk;
							A[ijk][j] += A[i][k]*t;	
						}
					}
				}
			}	
		}		
	}
	k = n;
	if (A[md][n-1] == 0.0) {
		ier = k;
		ip[n-1] = 0;
	}
	
	return (ier);

}

void solb(const int n, double **A, int ml, int mu, double *b, int *ip)
{
	
/*-----------------------------------------------------------------------
	Solution of linear system A*x = b

	Input:
		n 		order of matrix
		A		triangularized matrix obtained from decb
		ml		lower bandwidth of A (diagonal not counted)
		mu		upper bandwidth of A (diagonal not counted)
		b		right hand side vector
		ip		pivot vector obtained from dec
		
	Do not use if decb has set ier != 0
	
	Output:
		b		solution vector, x

-----------------------------------------------------------------------*/
	  
	int m, kb, nm1, md, md1, mdl, mdm, imd, kmd, lm;
	double t;

	md = ml + mu;
	md1 = md + 1;
	mdm = md - 1;
	nm1 = n - 1;	  
	if (n != 1) {
		if (ml != 0) {
			for (int k = 0; k < nm1; k++) {
				m = ip[k];
				t = b[m];
				b[m] = b[k];
				b[k] = t;
				mdl = min(ml, n - k - 1) + md;
				for (int i = md1; i <= mdl; i++) {
					imd = i + k - md;
					b[imd] += A[i][k]*t;
				}
			}
		}
		for (int k = 0; k < nm1; k++) {
			kb = n - k - 1; 
			b[kb] = b[kb]/A[md][kb];
			t = -b[kb];
			kmd = md - kb;
			lm = max(0, kmd); 
			for (int i = lm; i <= mdm; i++) {
				imd = i - kmd; 
				b[imd] += A[i][kb]*t;
			}
		}
	}
	b[0] = b[0]/A[md][0]; 
	
	return;

}

int decbc(const int n, double **AR, double **AI, int ml, int mu, int *ip)
{
	  
/*-----------------------------------------------------------------------
	Matrix Triangularization by Gaussian Elimination of a banded complex
		matrix with lower bandwidth ml and upper bandwidth mu
	
	Input:
		n			order of the original matrix A
		AR, AI		contains the matrix in band storage.
					The columns of the matrix are stored in the columns
					of AR (real part) and AI (imaginary part) and the 
					diagonals of the matrix are stored in rows ml through
					2*ml+mu of AR and AI
		ml			lower bandwidth of A (diagonal is not counted)
		mu			upper bandwidth of A (diagonal is not counted)
				
	Output:
		AR, AI		an upper triangular matrix in band storage and the 
						multipliers which were used to obtain it
		ip			index vector of pivot indices
		ip[n-1]		(-1)^(number of interchanges) or 0
		ier 		0 if matrix A is nonsingular, or k if found to be 
						singular at stage k
			
	Use solbc to obtain solution of linear system
	
	determ(A) = ip[n-1]*A[md][0]*A[md][1]*...*A[md][n-1] with md = ml+mu
	
	If ip[n-1] = 0, A is singular, solbc will divide by zero
	
	Reference:
		This is a modification of:
		C.B. Moler, Algorithm 423, Linear Equation Solver,
		C.A.C.M. 15 (1972), p. 274.
-----------------------------------------------------------------------*/	  

	int kp1, m, nm1, k, md, md1, mdl, mm, jk, ijk, ju, ier;
	double tr, ti, den, prodr, prodi;
	  
	ier = 0;
	ip[n-1] = 1;
	md = ml + mu;
	md1 = md + 1;
	ju = 0;
	if ((n != 1) && (ml != 0)) {
		if (n >= mu+2)
			for (int j = mu + 1; j < n; j++) 
				for (int i = 0; i < ml; i++) { 
					AR[i][j] = 0.0;
					AI[i][j] = 0.0;
				}
		nm1 = n - 1;
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;
			m = md;
			mdl = min(ml, n-k-1) + md;
			for (int i = md1; i <= mdl; i++) { 
				if ((fabs(AR[i][k]) + fabs(AI[i][k])) >
					(fabs(AR[m][k]) + fabs(AI[m][k]))) 
					m = i;
			}
			ip[k] = m + k - md;
			tr = AR[m][k];
			ti = AI[m][k];
			if (m != k) {
				ip[n-1] = -ip[n-1];
				AR[m][k] = AR[md][k];
				AI[m][k] = AI[md][k];
				AR[md][k] = tr;
				AI[md][k] = ti;
			}
			if ((fabs(tr) + fabs(ti)) == 0.0) {
				ier = k;
				ip[n-1] = 0;
				return (ier);
			}
			den = tr*tr + ti*ti;
			tr = tr/den;
			ti = -ti/den;
			for (int i = md1; i <= mdl; i++) { 
				prodr = AR[i][k]*tr - AI[i][k]*ti;
				prodi = AI[i][k]*tr + AR[i][k]*ti;
				AR[i][k] = -prodr;
				AI[i][k] = -prodi;
			}
			ju = min(max(ju, mu+ip[k]+1), n); 
			mm = md;
			if (ju >= kp1) {
				for (int j = kp1; j < ju; j++) {
					m--;
					mm--;
					tr = AR[m][j];
					ti = AI[m][j];
					if (m != mm) {
						AR[m][j] = AR[mm][j];
						AI[m][j] = AI[mm][j];
						AR[mm][j] = tr;
						AI[mm][j] = ti;
					}
					if ((fabs(tr) + fabs(ti)) == 0.0) {
					}
					else if (ti == 0.0) {
						jk = j - k;
						for (int i = md1; i <= mdl; i++) { 
							ijk = i - jk;
							prodr = AR[i][k]*tr;
							prodi = AI[i][k]*tr;
							AR[ijk][j] += prodr;
							AI[ijk][j] += prodi;
						}
					}
					else if (tr == 0.0) {
						jk = j - k;
						for (int i = md1; i <= mdl; i++) { 
							ijk = i - jk;
							prodr = -AI[i][k]*ti;
							prodi = AR[i][k]*ti;
							AR[ijk][j] += prodr;
							AI[ijk][j] += prodi;
						}
					}
					else {
						jk = j - k;
						for (int i = md1; i <= mdl; i++) { 
							ijk = i - jk;
							prodr = AR[i][k]*tr - AI[i][k]*ti;
							prodi = AI[i][k]*tr + AR[i][k]*ti;
							AR[ijk][j] += prodr;
							AI[ijk][j] += prodi;
						}
					}
				}
			}
		}	
	}	
	k = n;
	if ((fabs(AR[md][n-1]) + fabs(AI[md][n-1])) == 0.0) { 
		ier = k;
		ip[n-1] = 0;
	}

	return (ier);

} // decbc

void solbc(const int n, double **AR, double **AI, int ml, int mu,
	double *br, double *bi, int *ip)
{
	  
/*-----------------------------------------------------------------------
	Solution of linear system A*x = b
		version banded and complex
	
	Input:
		n			order of matrix
		AR, AI		triangularized matrix obtained from decb
						(real and imaginary parts)
		ml			lower bandwidth of A (diagonal is not counted)
		mu			upper bandwidth of A (diagonal is not counted)
		br, bi		right hand side vector (real and imaginary parts)
		ip			pivot vector obtained from decbc
		
	Do not use if decbc has set ier != 0
	
	Output:
		br, bi		solution vector, x (real and imaginary parts)

-----------------------------------------------------------------------*/
	  
	int nm1, m, kb, md, md1, mdl, mdm, imd, kmd, lm;
	double den, prodr, prodi, tr, ti;
  
	md = ml + mu;
	md1 = md + 1;
	mdm = md - 1;
	nm1 = n - 1;	  
	if (n != 1) {
		if (ml != 0) {
			for (int k = 0; k < nm1; k++) {
				m = ip[k];
				tr = br[m];
				ti = bi[m];
				br[m] = br[k];
				bi[m] = bi[k];
				br[k] = tr;
				bi[k] = ti;
				mdl = min(ml, n-k-1) + md;
				for (int i = md1; i <= mdl; i++) { 
					imd = i + k - md;
					prodr = AR[i][k]*tr - AI[i][k]*ti;
					prodi = AI[i][k]*tr + AR[i][k]*ti;
					br[imd] += prodr;
					bi[imd] += prodi;
				}
			} 
		}
 		for (int k = 0; k < nm1; k++) {
			kb = n - k - 1;
			den = AR[md][kb]*AR[md][kb] + AI[md][kb]*AI[md][kb];
			prodr = br[kb]*AR[md][kb] + bi[kb]*AI[md][kb];
			prodi = bi[kb]*AR[md][kb] - br[kb]*AI[md][kb];
			br[kb] = prodr/den;
			bi[kb] = prodi/den;
			tr = -br[kb];
			ti = -bi[kb];
			kmd = md - kb;
			lm = max(0, kmd);
			for (int i = lm; i <= mdm; i++) { 
				imd = i - kmd;
				prodr = AR[i][kb]*tr - AI[i][kb]*ti;
				prodi = AI[i][kb]*tr + AR[i][kb]*ti;
				br[imd] += prodr;
				bi[imd] += prodi;
			}
		}
		den = AR[md][0]*AR[md][0] + AI[md][0]*AI[md][0];
		prodr = br[0]*AR[md][0] + bi[0]*AI[md][0];
		prodi = bi[0]*AR[md][0] - br[0]*AI[md][0];
		br[0] = prodr/den;
		bi[0] = prodi/den;
	}
	
	return;

}


void elmhes(const int n, int low, int igh, double **A, int *inter)
{

/* --------------------------------------------------------------------
	This subroutine is a translation of the algol procedure elmhes,
	Num. Math. 12, 349-368(1968) by Martin and Wilkinson.
	Handbook for Auto. Comp., Vol.II-Linear Algebra, 339-358(1971).

	Given a real general matrix, this subroutine reduces a submatrix 
	situated in rows and columns low through igh to upper Hessenberg 
	form by stabilized elementary similarity transformations.

	Input:
		n 		order of the matrix;

		low, igh	integers determined by the balancing subroutine balanc. 
	  				If balanc has not been used, set low=0, igh=n;

		A 		the input matrix.

	Output:

		A		contains the Hessenberg matrix. The multipliers which 
	  			were used in the reduction are stored in the remaining 
				triangle under the Hessenberg matrix;

		inter	contains information on the rows and columns
				interchanged in the reduction. Only elements low 
				through igh are used.

	Questions and comments should be directed to B. S. Garbow,
	Applied Mathematics Division, Argonne National Laboratory

------------------------------------------------------------------*------*/

	int ii, la, kp1, mm1, mp1;
	double x, y;

	la = igh - 2;
	kp1 = low + 1;
	if (la < kp1) return;
	for (int m = kp1; m <= la; m++) {
		mm1 = m - 1;
		x = 0.0;
		ii = m;
		for (int j = m; j < igh; j++) {
			if (fabs(A[j][mm1]) > fabs(x)) {
				x = A[j][mm1];
				ii = j;
			}
		}
		inter[m] = ii;
		if (ii != m) {		
//    :::::::::: interchange rows and columns of a ::::::::::
			for (int j = mm1; j < n; j++) {
				y = A[ii][j];
				A[ii][j] = A[m][j];
				A[m][j] = y;
			}
			for (int j = 0; j < igh; j++) {
				y = A[j][ii];
				A[j][ii] = A[j][m];
				A[j][m] = y;
			}
		}  
//    :::::::::: end interchange ::::::::::
		if (x != 0.0) {
			mp1 = m + 1;
			for (int i = mp1; i < igh; i++) {
				y = A[i][mm1];
				if (y == 0.0) return;
				y = y/x;
				A[i][mm1] = y;
				for (int j = m; j < n; j++) A[i][j] -= y*A[m][j];
				for (int j = 0; j < igh; j++) A[i][m] += y*A[j][i];
			}
		}
	}

	return;

} //elmhes
