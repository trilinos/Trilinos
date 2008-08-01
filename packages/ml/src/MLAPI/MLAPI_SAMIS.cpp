/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#ifdef COMPILE_ME

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI.h"
#include "Epetra_BlockMap.h"
#include "Epetra_VbrMatrix.h"

namespace MLAPI {
namespace SAMIS {

// File reader written by M. Brezina; imported into MLAPI by M. Sala 
// on 07-Mar-05.

/* commented out on 14-Mar-05 to compile on atlantis (SGI64)
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
*/
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

using namespace std;

#define  KER_FILE        "ker.dat" 
#define  MTX_FILE        "mtx.dat" 

///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////

int 	prefetchMtx(
           const int  ascii, 
           int       &nf, 
           int       &m, 
           int       &n, 
           int       &s,
           const char myMtxFileName[]) 
{
       /**********************************************************************
	* Open the matrix data file and extract the matrix dimensions.
	* Close the file. 
	*
	* If myMtxFileName != NULL, fetch from the file of given name,
	* otherwise read from predefined filename.
	*
	* Arguments: 
	*   ascii         -in- 0=input file in binary format, 1=in ASCII form
	*   nf            -OU- dofs/node
	*   m             -OU- row dimension of matrix A 
	*   n             -OU- col dimension of matrix A 
	*   s             -OU- number of stored (block) entries in A 
	*   myMtxFileName -in- if != NULL, fetch from this file 
	*
	* Returns: 
	*     0  if all OK
	*     1  wrong binary format (endian) suspected 
	*     2  other severe error
        **********************************************************************/

	int      fd, ibuf[4], pad0, pad1; 
        size_t   szInt = sizeof(int);
        ssize_t  er;
	ifstream fdS; 
        char     fileName[128]; 
 
        if (myMtxFileName != NULL) { 
            sprintf(fileName, "%s", myMtxFileName);
        } else { 
            sprintf(fileName, "%s", MTX_FILE);
        }

	switch (ascii) { 
           case  0:/* binary file, should check endian type */
	            fd = open(fileName, O_RDONLY ); 
                    if (fd < 1) { 
                        cout   << "[SS] prefetchMtx: could not open file\n";
                        cout   << "[SS] prefetchMtx: could not open file\n";
                        return 2; 
                    }
//
//                  first read the record containing 4 scalar params nf,m,n,s
//
                    er = read(fd, &pad0, szInt); 
                    if (er != (ssize_t)szInt) {
                        cout   << "[SS] prefetchMtx: read pad failed\n";
                        return 2; 
                    }
                    if (pad0 != 4*(int)szInt) { 
                        cout   << "[EE] prefetchMtx: incorrect bin format\n";
                        return 1; 
                    }

		    er = read(fd, ibuf, 4*szInt);  

                    if (er != (ssize_t) (4*szInt)) {
                        cout   << "[SS] prefetchMtx: read failed\n";
                        cout         << "[SS] prefetchMtx: read failed\n"; 
                        return 2; 
                    }
                    nf = ibuf[0]; 
                    m  = ibuf[1]; 
                    n  = ibuf[2]; 
                    s  = ibuf[3]; 
                    
                    er = read(fd, &pad1, szInt); 

                    if (er != (ssize_t)szInt) {
                        cout   << "[SS] prefetchMtx: read pad failed\n";
                        return 2; 
                    }
                    if (pad1 != pad0) { 
                        cout   << "[SS] prefetchMtx: file corrupt ?\n";
                        return 2; 
                    }
                    
                    close(fd); 
		    break; 
           case  1:/* ASCII file in the (nf, m, n, s; x; r; v) format  */
                    fdS.open(fileName); 
                    if (fdS.fail()) {
                        cout   << "[SS] prefetchMtx: could not open file " 
                               << fileName << "\n";
                        return 2; 
                    }
//
//                  ascii format starts with 4 lines containing the 
//                  4 scalar params nf,m,n,s
//
                    fdS >> nf >> m >> n >> s; 
	            fdS.close(); 
		    break; 
	   default: /* invalid entry */
                    cout   << "[SS] prefetchMtx: invalid ascii\n";  
	            return 2; 
	}

  	return 0; 
}

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

int     fetchMtx(
           const int  ascii, 
           const int  nf, 
           const int  m, 
           const int  n, 
           const int  s, 
           int       *x, 
           int       *r, 
           double    *v,
           const char myMtxFileName[]
        ) 
{
       /**********************************************************************
	* Open the matrix data file and extract the matrix.
	* The dimensions are assumed known and the matrix itself already 
	* allocated before calling this function.
	*
	* If myMtxFileName != NULL, fetch from the file of given name,
	* otherwise read from predefined filename.
	* 
	* To extract the dimensions, use prefetchMtx(). 
	*
	* Arguments: 
	*   ascii     -in- 0=input file in binary format, 
	*                  1=in ASCII format (nf,m,n,s; x; r;v)
	*                  2=in ASCII format Aij format (for scalar mt only)
	*   nf        -in- dofs/node
	*   m         -in- row dimension of matrix A 
	*   n         -in- col dimension of matrix A 
	*   s         -in- number of stored (block) entries in A 
	*   x[n+1]    -OU- pointers to row info for each column (F77 indexing)
	*   r[s]      -OU- row indices                          (F77 indexing) 
	*   v[nf*nf*s]-OU- (block) entries                      (F77 indexing) 
	*   myMtxFileName -in- if != NULL, fetch from this file
        *
	* Returns: 
	*     0  if all OK
	*     1  wrong binary format (endian) suspected 
	*     2  other severe error
        **********************************************************************/

	int     fd, ibuf[4], i; 
        size_t  szInt = sizeof(int);
        size_t  szDbl = sizeof(double);
        size_t  sz;
        ssize_t er;
	int     pad0, pad1; 
	int     recordMark0[]=
                       {4*(int)sizeof(int),-1, -2, -3, -4, 4*(int)sizeof(int)};
	int     recordMark1[6]; // 6 ints (pad + 4 + pad) storing the file mark
        char    fileName[128]; 
        ifstream fdS;
 
        if (myMtxFileName != NULL) { 
            sprintf(fileName, "%s", myMtxFileName);
        } else { 
            sprintf(fileName, "%s", MTX_FILE);
        }

        if (nf < 1   || m < 1 || n < 1 || s < 1) { 
            cout << "[SS] fetchMtx: invalid mtx params: nf=" << nf
                 << " m=" << m << " n=" << n << " s=" << s << endl;
            abort();
        }

	switch (ascii) { 
           case  0:/* binary file, should check endian type */
	            fd = open(fileName, O_RDONLY ); 
                    if (fd < 1) { 
                        cout   << "[SS] fetchMtx: could not open file\n";  
                        return 2; 
                    }
//
//                  first read the record containing 4 scalar params nf,m,n,s
//
                    sz = 4 * szInt; 
                    er = read(fd, &pad0, szInt); 

                    if (er != (ssize_t)szInt) {
                        cout   << "[SS] fetchMtx: read pad failed\n";
                        return 2; 
                    }
                    if (pad0 != (int)sz) { 
                        cout   << "[SS] fetchMtx: mtx in wrong bin format\n";
                        close(fd); 
                        return 1; 
                    }

		    er = read(fd, ibuf, sz);  

                    if (er != (ssize_t)sz) {
                        cout   << "[SS] fetchMtx: read failed\n"; 
                        return 2; 
                    }
	            if (ibuf[0]!=nf || ibuf[1]!=m || ibuf[2]!=n || ibuf[3]!=s) {
                        cout   << "[SS] fetchMtx: incompatible mtx file !\n";
                        return 2; 
                    }
                    
                    er = read(fd, &pad1, szInt); 

                    if (er != (ssize_t)szInt) {
                        cout   << "[SS] fetchMtx: read pad1 failed\n";
                        return 2; 
                    }
                    if (pad1 != pad0) { 
                        cout   << "[SS] fetchMtx: invalid pad1\n"; 
                        return 2; 
                    }
//
//                  check filemark (read both pads here at he same time)
//
                    er = read(fd, recordMark1, 6*szInt); 

                    for (i=0; i<6; i++) { 
                         if (recordMark1[i] != recordMark0[i]) { 
                             cout   << "[SS] fetchMtx: invalid filemark!\n";
                             return 2; 
                         }
                    }
//
//                  read the x vector 
//
                    er = read(fd, &pad0, szInt      ); 
                    er = read(fd, x,     (n+1)*szInt); 
                    er = read(fd, &pad1, szInt      ); 
			
                    if (s != x[n] - x[0]) { 
                        cout   << "[EE] fetchMtx: param s != x[n] - x[0]\n";
                        return 2; 
                    }
//
//                  check filemark (read both pads here at he same time)
//
                    er = read(fd, recordMark1, 6*szInt); 
                    for (i=0; i<6; i++) { 
                         if (recordMark1[i] != recordMark0[i]) { 
                             cout   << "[SS] fetchMtx: invalid filemark!\n";
                             return 2; 
                         }
                    }
//
//                  read the r vector: 
//
	            er = read(fd, &pad0, szInt  );
                    er = read(fd, r,     s*szInt); 
	            er = read(fd, &pad1, szInt  );
//
//                  check filemark (read both pads here at he same time)
//
                    er = read(fd, recordMark1, 6*szInt); 
                    for (i=0; i<6; i++) { 
                         if (recordMark1[i] != recordMark0[i]) { 
                             cout   << "[SS] fetchMtx: invalid filemark!\n";
                             return 2; 
                         }
                    }
//
//                  read the v vector: 
//
	            er = read(fd, &pad0, szInt          );
                    er = read(fd, v,     (nf*nf*s)*szDbl); 
	            er = read(fd, &pad1, szInt          );

                    close(fd); 
		    break; 
           case  1:/* ASCII file in the (nf, m, n, s; x; r; v) format */
                    fdS.open(fileName); 
                    if (fdS.fail()) { 
                        cout << "[SS] fetchMtx: could not open file\n";
                        return 2;
                    }
                    fdS >> ibuf[0]; 
                    fdS >> ibuf[1]; 
                    fdS >> ibuf[2]; 
                    fdS >> ibuf[3]; 
	            if (ibuf[0]!=nf || ibuf[1]!=m || ibuf[2]!=n || ibuf[3]!=s) {
                        cout   << "[SS] fetchMtx: incompatible mtx file !\n";
                        return 2; 
                    }
                    for(i=0; i<n+1;     i++) fdS >> x[i];
                    for(i=0; i<s;       i++) fdS >> r[i];
                    for(i=0; i<s*nf*nf; i++) fdS >> v[i];

                    fdS.close(); 
		    break; 
	   default:/* invalid entry */
                    cout   << "[SS] fetchMtx: invalid ascii\n";  
	            return 2; 
	}

        cout << "[II] extracted matrix from file " << fileName << endl;
	
  	return 0; 
}

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

int 	prefetchKers(
           const int  ascii, 
           int       &nDof, 
           int       &nNod, 
           int       &nKer,
           const char myKerFileName[]
        )
{
       /**********************************************************************
	* Open the kernel data file and extract the dimensions.
	* Close the file. 
	*
	* If myKerFileName is not NULL, use that filename instead of std.
	*
	* Arguments: 
	*   ascii         -in- 0=input file in binary format, 
	*                      1=input in ASCII form
	*                      2=input in ASCII form (for consistency with mtx)
	*   nNod          -OU- number of nodes 
	*   nDof          -OU- number of dofs/node
	*   nKer          -OU- number of kernel vectors 
	*   myKerFileName -in- filename to prefetch from instead of predef.
        *
	* Returns: 
	*     0  if all OK
	*     1  wrong binary format (endian) suspected 
	*     2  other severe error
        **********************************************************************/

	int      fd; 
        int      pad0, pad1; 
        ssize_t  er; 
        size_t   szInt=sizeof(int); 
        char     fileName[128]; 
        ifstream fdS;

        if (myKerFileName != NULL) { 
            sprintf(fileName, "%s", myKerFileName);
        } else { 
            sprintf(fileName, "%s", KER_FILE);
        }

	switch (ascii) {
           case  0:/* Kernels are stored in a binary file */ 
                    fd = open(fileName, O_RDONLY); 
                    if (fd <1) { 
                        cout   << "[SS] prefetchKers: could not open file\n";
                        return 2; 
                    }
//
//                  fetch the scalar data from file
//
                    er=read(fd, &pad0, szInt); 
                    er=read(fd, &nNod, szInt); 
                    er=read(fd, &pad1, szInt); 

                    if (er != (ssize_t)szInt) {
                        cout << "[SS] prefetchKers: reading failed\n";
                        abort();
                    }

                    if (pad0 != (int)szInt) {
                        cout   << "[SS] prefetchKers: wrong bin format\n";
                        close(fd); 
                        return 1; 
                    }

                    er=read(fd, &pad0, szInt); 
                    er=read(fd, &nDof, szInt); 
                    er=read(fd, &pad1, szInt); 

                    er=read(fd, &pad0, szInt); 
                    er=read(fd, &nKer, szInt); 
                    er=read(fd, &pad1, szInt); 

                    close(fd); 
                    break; 
           case  1:/* Kernels are stored in a ASCII file */ 
           case  2:/* Kernels are stored in a ASCII file */ 
                    fdS.open(fileName);
                    fdS >> nNod >> nDof >> nKer; 
                    fdS.close();
                    break; 
           default:/* unsupported value of ascii */
                    cout   << "[SS] prefetchKers: invalid ascii="
                             << ascii << endl;
                    exit(1);
	}

        cout << "[DD]  G0 prefetchKer: got nNod=" << nNod
             << ", nDof=" << nDof << ", nKer=" << nKer << endl;

	if (nKer <1 && nNod < 1 || nDof <1) { 
            cout   << "[SS] prefetchKers: read invalid dimensions\n"; 
            cout   << " ...    nNod="<<nNod<<" nDof="<<nDof
                   << " nKer="<<nKer<<endl; 
            return 2; 
        }

	return 0; 
}

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

int	fetchKers(
        const int  ascii,          /* -in- file format 0=bin 1,2= ascii      */
        const int  nf,             /* -in- number of dof/node                */
        const int  nNod,           /* -in- number of nodes                   */
        const int  nKer,           /* -in- number of kernel vectors          */
        const int  limKer,         /* -in- if >= 0, limit #ker actually used */
        double    *xKer,           /* -ou- ker. vectors in one flat array    */
        const char myKerFileName[] /* -in- filename to read from or NULL(std)*/
        )
{
       /**********************************************************************
        * Fetch the kernel vectors from the file defined by KER_FILE.
        * This reads the block-version of the kernels; i.e., 
        * the kernel vectors may have been padded with zero values 
        * for the artificial degrees of freedom to ensure that 
        * the number of dof pre node is constant. 
        *
        * If myKerFileName is not NULL, use that filename instead of std.
        *
        * The file is assumed in F77 unformatted format, so we have to 
        * take care of handling the record pads here.
        * The file consits of the following F77 records: 
        *   1.    nnods
        *   2.    nnf 
        *   3.    nKers
        *         do i=1, nKers
        *   4,...    kernel_vector[i]
        *         enddo
        *
        * Arguments: 
        *   ascii         -in- 0=input file in binary format, 
        *                      1=input in ASCII form
        *                      2=input in ASCII form (for consistency with mtx)
        *   nNod          -in- number of nodes       (from prefetchKers())
        *   nf            -in- number dof/node       (from prefetchKers())
        *   nKer          -in- number kernel vectors (from prefetchKers())
        *   xKer[]        -OU- ker vectors in one flat array allocated outside
        *   myKerFileName -in- filename to prefetch from instead of predef.
        *
        * Returns: 
        *     0  if all OK
        *     1  wrong binary format (endian) suspected 
        *     2  other severe error
        **********************************************************************/
        int      i, pad0, pad1, f_ker, nnf, nnods, nkers; 
        size_t   szInt  = sizeof(int); 
        size_t   szDbl  = sizeof(double); 
        size_t   szKer; 
        ssize_t  rs; 
        char     fileName[128]; 
        ifstream fdS;

        if (myKerFileName != NULL) {
            sprintf(fileName, "%s", myKerFileName);
        } else {
            sprintf(fileName, "%s", KER_FILE);
        }

        switch(ascii) {
           case  0:/* Kernels are stored in a BINARY file */ 
                    f_ker = open(fileName, O_RDONLY); 
            
            	    if (f_ker < 0) { 
                        cout   << "[SS] fetchKer: could not open file "
                               << fileName <<endl;
                        return 2; 
                    }
             
                    rs=read(f_ker, &pad0,   szInt   ); 
                    rs=read(f_ker, &nnods,  szInt   );  
                    rs=read(f_ker, &pad1,   szInt   ); 
            
                    if (pad0 != pad1 ) { 
                        cout   << "[SS] fetchKer: file corrupt\n";
                        return 2;  
                    }
            
                    if ( pad0 !=(int)szInt ) { 
                        cout   << "[EE] fetchKer: wrong-endian format ?\n";
                        return 1;  
                    }
                    if (nnods != nNod) { 
                        cout   << "[SS] fetchKer: read nnods=" << nnods
                               << " incosistent W/" << " nNod=" << nNod << endl;
                        return 2;  
                    }
            
                    rs=read(f_ker, &pad0, szInt   ); 
                    rs=read(f_ker, &nnf,  szInt   );  
                    rs=read(f_ker, &pad1, szInt   ); 
            
                    if (nnf != nf) {
                        cout   << "[SS] fetchKer: read nnf="
                               << nnf << " incosistent\n";
                        return 2;  
                    }
                    if (pad0 != pad1 ) { 
                        cout   << "[SS] fetchKer: file corrupt\n";
                        return 2;  
                    }
            
                    rs=read(f_ker, &pad0,   szInt); 
                    rs=read(f_ker, &nkers,  szInt);  
                    rs=read(f_ker, &pad1,   szInt); 
            
                    if (nkers != nKer) {
                        if (nKer == limKer) {
                            cout   << "[WW] fetchKer: limiting nkers="
                                   << nkers<< " to " << nKer << endl;
                        } else { 
                            cout   << "[SS] fetchKer: read nkers="
                                   << nkers<<" incosistent\n";
                            return 2;  
                        }
                    }
                    if (pad0 != pad1 ) {
                        cout   << "[SS] fetchKer: file corrupt\n";
                        return 2;  
                    }
            
                    szKer = nf*nNod*szDbl; 
            
                    for (i=0; i<nKer; i++) { 
                         
                         rs=read(f_ker, &pad0, szInt); 
            
                         if ( pad0 != (int)szKer ) {  
                              cout   << "[SS] fetchKer: file corrupt\n";
                              return 2;  
                         }
            
                         rs=read(f_ker, &xKer[i*nf*nNod],  pad0);  
                         rs=read(f_ker, &pad1,             szInt); 
                    }
                    if (pad0 != pad1 ) { 
                        cout   << "[SS] fetchKer: file corrupt\n";
                        return 2;  
                    }
            
                    close(f_ker); 
                    break;
           case  1:/* Kernels are stored in a ASCII file */ 
           case  2:/* Kernels are stored in a ASCII file */ 
                    fdS.open(fileName);
                    if(fdS.fail()) { 
                        cout << "[SS] fetchKer: failed to open file\n";
                        abort();
                    }

                    fdS >> nnods >> nnf >> nkers;

                    if (nnods!=nNod || nnf != nf ) {
                        cout << "[SS] fetchKer: header inconsistency\n";
                        abort();
                    }
                    if (nkers != nKer) { 
                        if (nKer == limKer) {
                            cout   << "[WW] fetchKer: limiting nkers="
                                   << nkers<< " to " << nKer << endl;
                        } else { 
                            cout   << "[SS] fetchKer: read nkers="
                                   << nkers<<" incosistent\n";
                            return 2;  
                        }
                    }

                    for (i=0; i<nKer; i++) {
                         for (int j=0; j <nNod*nf; j++) { 
                              fdS >> xKer[i*nf*nNod + j]; 
                         }
                    }
                    fdS.close();
                    break;
           default:
                    cout << "[SS] fetchKer: illegal ascii=" << ascii << endl;
                    abort();
        }

        return 0; 
}

} // namespace SAMIS

// ====================================================================== 
void ReadSAMISMatrix(const char *filen, Operator& A, int& NumPDEEqns) 
{
  if (GetNumProcs() != 1)
    ML_THROW("SAMIS interface works in serial only!", -1);

  /**************************************************************************
   * Example of usage: fetches the matrix from SAMISdat(AMG) binary 
   * file (block version) from a file specified by arg. filen.
   * If filen==NULL, use default name sepcified by MTX_FILE.
   **************************************************************************/

  int     nf, m, n, s;
  int     ascii = 0; // want binary
  int     er, blkSz;

  er = MLAPI::SAMIS::prefetchMtx(ascii, nf, m, n, s, filen);

  if (er)
    ML_THROW("prefetchMtx return value = " + GetString(er), -1);

  if (GetPrintLevel())
    fprintf(stderr,"[DD] read bin mtx: nf=%8d m=%8d n=%8d s=%9d\n", 
            nf, m, n, s);

  // build up a block matrix with n blocks, each with nf dof's
  // on one processor only.

  Epetra_BlockMap BlockMap(n, nf, 0, GetEpetra_Comm());
  Epetra_VbrMatrix* VbrMatrix = new Epetra_VbrMatrix(Copy, BlockMap, 0);

  blkSz = nf*nf;
  vector<int>    x(n + 1);
  vector<int>    r(s);
  vector<double> v(s * blkSz);

  er = MLAPI::SAMIS::fetchMtx(ascii, nf, m, n, s, &x[0], &r[0], &v[0], filen);

  if (er)
    ML_THROW("fetchMtx return value = " + GetString(er), -1);
  
  if (x[0] != 1) {
    fprintf(stderr,"[SS] input assumed from F77, so x[0]=%d impossible\n",
            x[0]);
    abort();
  }

  vector<double> Values(blkSz);

  for (int i = 0; i < n ; i++) {
    
    int NumEntries = x[i + 1] - x[i];
    vector<int> Indices(NumEntries);
    int count = 0;
    for (int j=x[i]-1; j < x[i+1]-1; j++)
      Indices[count++] = r[j] - 1;

    VbrMatrix->BeginInsertGlobalValues(i, NumEntries, &Indices[0]);

    for (int j = x[i] - 1; j < x[i+1] - 1; j++) {
      for (int k = 0 ; k < blkSz ; ++k )
        Values[k] = v[blkSz * j + k];

      VbrMatrix->SubmitBlockEntry(&Values[0],nf,nf,nf);
    }

    VbrMatrix->EndSubmitEntries();
  }

  VbrMatrix->FillComplete();

  Space MatrixSpace(nf * n);

  A.Reshape(MatrixSpace, MatrixSpace, VbrMatrix, true);

  NumPDEEqns = nf;

} // ReadSAMISMatrix

// ====================================================================== 
void ReadSAMISKernel(const char *myKerFileName, MultiVector& A, 
         const int limKer)
{

  if (GetNumProcs() != 1)
    ML_THROW("SAMIS interface works in serial only!", -1);

  int ascii = 0;
  int nDof, nNod, nKer;
  int err;

  err = MLAPI::SAMIS::prefetchKers(ascii, nDof, nNod, nKer, myKerFileName);
  if (err)
    ML_THROW("prefetchKers return value = " + GetString(err), -1);

  vector<double> data(nDof * nNod * nKer);

  err = MLAPI::SAMIS::fetchKers(ascii, nDof, nNod, nKer, limKer, 
                                &data[0], myKerFileName);
  if (err)
    ML_THROW("fetchKers return value = " + GetString(err), -1);

  Space VectorSpace(nDof * nNod);

  A.Reshape(VectorSpace, nKer);
  for (int v = 0 ; v < nKer ; ++v)
    for (int i = 0 ; i < nDof * nNod ; ++i)
      A(i, v) = data[i + v * nDof * nNod];

  return;
  
} // ReadSAMISKernel

} // namespace MLAPI

#endif // HAVE_ML_MLAPI
#endif
