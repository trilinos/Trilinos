/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/****************************************************************************
 * PETSc_LinSysCore.C              Mark F. Adams                     May 2000
 *
 *   LinearSystemCore derived class for PETSc
 *  
 * Limitations:
 *  1) Assumes symmetric structure of matrix. 
 *  2) Assumes SYMMETRIC MATRICES to apply BCs (easy to fix).
 *      - need a good PETSc getcol() to fix (enforceEssentialBC)
 *      - need a good PETSc zerocol() or getcol() to zero columns
 *      - slow getcol (MatGetColumnVector) is available in PETSc but
 *        2.0.28 has a bug, Barry Smith can give you a patch for 
 *        this.  Need to use it instead of MatGetRow() in one place.
 *
 * Notes,
 *  1) Parameters can be given in argument pairs "-option [value]", see the
 *     PETSc documentation for parameter specification (command line arguments)
 *  2) Constraints are handled with a simple "Uzawa" or Augmented Lagrange 
 *     solver that utilizes PETSc's solvers to solve the (normalized) 
 *     stiffness matrix part of the system.  Thus, if your operator is SPD 
 *     then you can use CG even though the total linear system is indefinite 
 *     and the (unnormalized) stiffness matrix part is only symmetric positive 
 *     semi-definite.
 *
 ****************************************************************************/
 
#include <math.h>
#include <stdlib.h>
#include <fei_iostream.hpp>
#include <assert.h>
#include <fei_Data.hpp>
#include <fei_Lookup.hpp>
#include <fei_LinearSystemCore.hpp>

#include <fei_PETSc_LinSysCore.hpp>

#ifdef FEI_USE_C_PETSC
extern "C" {
#endif

#include <petsc.h>
#include <petscksp.h>

#ifdef FEI_USE_C_PETSC
}
#endif

//petsc.h defines PetscMalloc(a,b) and PetscFree(a) to be macros that call
//other functions and pass those other functions the string macros __FUNCT__,
//__FILE__ and __SDIR__. Those string macros are string-literals but the
//function being called is expecting 'char*' arguments (i.e., non-const
//strings). Passing a string-literal to a function that's expecting a char*
//causes a warning on the Solaris C++ compiler. To get around that, we will
//re-define the PetscMalloc and PetscFree macros here to call the
//same underlying function, but pass strings which aren't string literals.

static char petsc_linsyscore_funct[] = "unknown";
static char petsc_linsyscore_file[] = "PETSc_LinSysCore.C";
static char petsc_linsyscore_sdir[] = "unknown";

#ifdef FEI_USE_C_PETSC
extern "C" {
#endif

#undef PetscMalloc
#define PetscMalloc(a,b)  (*PetscTrMalloc)((a),__LINE__, petsc_linsyscore_funct,petsc_linsyscore_file, petsc_linsyscore_sdir,(void**)(b))

#undef PetscFree
#define PetscFree(a)  (*PetscTrFree)((a),__LINE__, petsc_linsyscore_funct,petsc_linsyscore_file, petsc_linsyscore_sdir)

#undef CHKERRQ
#define CHKERRQ(n) if (n) {return PetscError(__LINE__,petsc_linsyscore_funct,petsc_linsyscore_file,petsc_linsyscore_sdir,n,0," ");}

#undef SETERRQ
#define SETERRQ(n,s) {return PetscError(__LINE__,petsc_linsyscore_funct,petsc_linsyscore_file,petsc_linsyscore_sdir,n,1,s);}

#undef SETERRQ2
#define SETERRQ2(n,s,a1,a2) {return PetscError(__LINE__,petsc_linsyscore_funct,petsc_linsyscore_file,petsc_linsyscore_sdir,n,1,s,a1,a2);}

#ifdef FEI_USE_C_PETSC
}
#endif

// PromCRVec ================================================================
class PromCRVec
{ 
 public:
  PromCRVec( Vec x = NULL, Vec p = NULL){ x_ = x; p_ = p; };
  ~PromCRVec(){ if( x_ ) VecDestroy(x_); if( p_ ) VecDestroy( p_ ); }
  //
  int Norm( NormType type, PetscReal *val );
  int Set( const double alpha );
  int AYPX( const double *alpha, PromCRVec *x );
  int AXPY( const double *alpha, PromCRVec *x );
  // data
  Vec x_;
  Vec p_; 
};

// PromCRMat ================================================================
class PromCRMat
{
 public:
  PromCRMat( Mat A = NULL, Mat C = NULL ){ A_ = A; C_ = C; }
  ~PromCRMat(){ if( A_ ) MatDestroy( A_ ); if( C_ ) MatDestroy( C_ ); }
  //
  int Mult( PromCRVec *x, PromCRVec *y );
  int ZeroEntries();
  // data
  Mat A_;
  Mat C_; 
};


// PETSc_OthBCData =========================================================
class PETSc_ZeroEquations
{
public:
  PETSc_ZeroEquations( int* localrows, int* sizes, int** idsarr, int len ) {
    len_ = len; 
    localrows_ = localrows;  sizes_ = sizes;
    idsarr_ = idsarr;
    next = NULL;
  }
  ~PETSc_ZeroEquations(){}
  int len_;
  int *sizes_;
  int *localrows_; 
  int **idsarr_;
  PETSc_ZeroEquations *next;
};

int printArr( FILE *file, const int *arr, const int n,
	      const char *header, const int step = 100 )
{
  fprintf(  file, header, n );      
  for ( int j = 0 ; j < n ; j++ ) {
    if ( (step > 0 && !(j%step)) || (!step && !j) ) {
      fprintf(  file, "\n");        
    }
    fprintf( file," %d", arr[j] );
  }
  
  fprintf( file, "\n");            
  return 0;
}

int printArr( FILE *file, const int *arr, const int n,
	      int extrad, const char *header,
	      const int step = 100 )
{
  fprintf(  file, header, n, extrad );      
  for ( int j = 0 ; j < n ; j++ ) {
    if ( (step > 0 && !(j%step)) || (!step && !j) ) {
      fprintf(  file, "\n");        
    }
    fprintf( file," %d", arr[j] );
  }
  
  fprintf( file, "\n");            
  return 0;
}

// PETSc_EssBCData =========================================================
class PETSc_EssBCData
{
public:
  PETSc_EssBCData(int* globalEqn, double* alpha, double* gamma, int len) {
    len_ = len; int ii,ierr;
    ierr = PetscMalloc(len*sizeof(int),&globalEqn_);  assert(ierr==0);
    ierr = PetscMalloc(len*sizeof(double),&alpha_);  assert(ierr==0);
    ierr = PetscMalloc(len*sizeof(double),&gamma_);
    if (ierr != 0) abort();
    for(ii=0;ii<len;ii++){
      globalEqn_[ii]=globalEqn[ii]; alpha_[ii]=alpha[ii]; gamma_[ii]=gamma[ii];
    }
    next = NULL;
  }
  ~PETSc_EssBCData(){
    PetscFree(globalEqn_); PetscFree(alpha_); PetscFree(gamma_);
  }
  int len_;
  int* globalEqn_;
  double* alpha_;
  double* gamma_;
  PETSc_EssBCData *next;
};

// PETSc_OthBCData =========================================================
class PETSc_OthBCData
{
public:
  PETSc_OthBCData( int* globalEqn, double* alpha, double* beta, double* gamma,
                   int len) {
    len_ = len; int ii,ierr;
    ierr = PetscMalloc(len*sizeof(int),&globalEqn_); assert(ierr==0);
    ierr = PetscMalloc(len*sizeof(double),&alpha_); assert(ierr==0);
    ierr = PetscMalloc(len*sizeof(double),&beta_); assert(ierr==0);
    ierr = PetscMalloc(len*sizeof(double),&gamma_);
    if (ierr != 0) abort();
    for(ii=0;ii<len;ii++){
      globalEqn_[ii] = globalEqn[ii]; alpha_[ii] = alpha[ii];
      beta_[ii] = beta[ii]; gamma_[ii] = gamma[ii];
    }
    next = NULL;
  }
  ~PETSc_OthBCData(){
    PetscFree(globalEqn_); PetscFree(alpha_);
    PetscFree(gamma_); PetscFree(beta_);
  }
  int len_;
  int* globalEqn_;
  double* alpha_;
  double* beta_;
  double* gamma_;
  PETSc_OthBCData *next;
};

#define ERR_RETURN(str){std::cerr<<str<<std::endl; return(-1);} //SETERRQ(1,str);}
#define PRINT if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d\n",thisProc_,__FUNC__,__LINE__)
#define PRINT1(a) if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d: %d\n",thisProc_,__FUNC__,__LINE__,a)
#define PRINTD(a) if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d: %e\n",thisProc_,__FUNC__,__LINE__,a)
#define PRINTS(a) if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d: %s\n",thisProc_,__FUNC__,__LINE__,a)
#define PRINT2(a,b) if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d: %d %d\n",thisProc_,__FUNC__,__LINE__,a,b)
#define PRINT3(a,b,c) if(verbose_>0)PetscPrintf((verbose_==1)?FEI_PETSc_Comm_:MPI_COMM_SELF,"[%d]%s line %d: %d %d %d\n",thisProc_,__FUNC__,__LINE__,a,b,c)

#undef CHKERRQ
#define CHKERRQ(e) if(e){fprintf(stdout,"[?]%s ERROR: line %d\n",__FUNC__,__LINE__); abort(); }

#undef TRUE
#undef FALSE
#define TRUE true
#define FALSE false

static char vector_name[] = "PETSc_Vector";
static char matrix_name[] = "PETSc_Matrix";
static double one = 1.0, zero = 0.0, mone = -1.0;
static char petscrc_name[] = "./.petscrc";

/***********************************************************************
 *                        PETSc_LinSysCore
 **********************************************************************/
int PETSc_LinSysCore::nActiveLSC_ = 0; 
int PETSc_LinSysCore::petsc_init_called_ = 0; 
// PETSc_LinSysCore::PETSc_LinSysCore() ======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::PETSc_LinSysCore"
PETSc_LinSysCore::PETSc_LinSysCore(MPI_Comm comm)
  : masterProc_(0),
    ndf_(0),
    dirty_guess_(FALSE),
    dirty_system_(TRUE),
    K_(NULL),
    x_(NULL),
    b_(NULL),
    init_guess_(NULL),
    Pen_stiff_void_(NULL),
    currentRHS_(-1),
    proc_gnode_(NULL),
    proc_globalEq_(NULL),
    localStartRow_(0),
    numLocalRows_(0),
    localEndRow_(0),
    numGlobalRows_(0),
    proc_lmEq_(NULL),
    localStartRowLm_(0),
    numLocalRowsLm_(0),
    localEndRowLm_(0),
    numGlobalRowsLm_(0),
    maxNodesPerCR_(0),
    proc_primEq_(NULL),
    localStartRowPrim_(0),
    numLocalRowsPrim_(0),
    localEndRowPrim_(0),
    numGlobalRowsPrim_(0),
    essbcdata_(NULL),
    othbcdata_(NULL),
    zeroEqs_(NULL),
    verbose_(0)
{
  static char dummy[] = "PETSc_LinSysCore";
  char* dummyptr = (char*)dummy;
  int argc = 1;  char **argv = &dummyptr; 
  FEI_PETSc_Comm_ = comm;
  MPI_Comm_size(comm, &numProcs_);
  MPI_Comm_rank(comm, &thisProc_);
  PRINT2(nActiveLSC_,petsc_init_called_);
  if( PetscInitializeCalled == PETSC_FALSE )  {
    PetscSetCommWorld( comm ); 
    PetscInitialize( &argc, &argv, petscrc_name, 0 );
  }

  numRHSs_ = 1;
#ifndef NDEBUG
  int ierr =
#endif
  PetscMalloc( sizeof(int)*numRHSs_,&rhsIDs_); assert(ierr==0);
  rhsIDs_[0] = 0;

  nActiveLSC_++;
  petsc_init_called_++;
}

//========DESTRUCTOR==========================================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::~PETSc_LinSysCore"
PETSc_LinSysCore::~PETSc_LinSysCore() 
{
  PRINT1(nActiveLSC_);

  if( K_ != NULL ) { delete K_; delete x_; delete init_guess_; }
  if( b_ != NULL ) {
    for( int ii = 0; ii < numRHSs_; ii++ ) if( b_[ii] != NULL ) delete b_[ii];
    PetscFree( b_ ); 
  }
  if( proc_gnode_ != NULL ) PetscFree( proc_gnode_ );
  if( proc_globalEq_ != NULL ) PetscFree( proc_globalEq_ );
  if( proc_lmEq_ != NULL ) PetscFree( proc_lmEq_ );
  if( proc_primEq_ != NULL ) PetscFree( proc_primEq_ );
  if( rhsIDs_ != NULL ) PetscFree( rhsIDs_ ); 
  
  nActiveLSC_--; 
//  petsc_init_called_--;
  // not clear when to call PetscFinalize!!!
//  if( nActiveLSC_ == 0 ) PetscFinalize();
}

//============================================================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::clone"
LinearSystemCore* PETSc_LinSysCore::clone() 
{
  PRINT;
  return( new PETSc_LinSysCore(FEI_PETSc_Comm_) );
}

// PETSc_LinSysCore::parameters ==============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::parameters"
int PETSc_LinSysCore::parameters( int numParams, const char*const* params ) 
{
  PRINT1(numParams);
  int ierr, ii, jj, done, i1, i2;
  char *iname, *value, buff1[128], buff2[128];

  for( ii = 0 ; ii < numParams ; ii++ ) {
    for( jj = i1 = i2 = done = 0, iname = value = NULL ; 
	 params[ii][jj] != '\0' ; 
	 jj++ ){
      if( !iname && params[ii][jj] == '-' ){ 
	iname = buff1; buff1[i1++] = params[ii][jj]; // start arg1
      }
      else if( params[ii][jj] == ' ' ) {
	if( iname && done == 0 ){ done++; buff1[i1++] = '\0';} // cap off arg1
	else if( value && done == 1 ){ 
	  done++; buff2[i2++] = '\0'; break; // cap off arg2
	}
	// else skip spaces
      }
      else if( iname && done == 0 ) buff1[i1++] = params[ii][jj]; // copy arg1
      else if( done == 1 ){ 
	if( value == NULL ) value = buff2; // start arg2
	buff2[i2++] = params[ii][jj]; // copy arg2
      }
    }
    if( iname && done == 0 ){ done++; buff1[i1++] = '\0';} // cap off arg1
    if( value && done == 1 ){ done++; buff2[i2++] = '\0';} // cap off arg2
    
    // call Petsc's OptionsSetValue
    if( done > 0 ) {
      if(verbose_)PetscPrintf( (verbose_==1) ? FEI_PETSc_Comm_ : MPI_COMM_SELF,
			       "\t[%d]%s call OptionsSetValue(%s,%s)\n",
			       thisProc_,__FUNC__, iname, value ? value : ""); 
      ierr = PetscOptionsSetValue( iname, value );  CHKERRQ(ierr);
    }
  }
  
  if( numParams == -1 ) {
    ierr = PetscOptionsSetValue( "-ksp_type", "gmres" );  CHKERRQ(ierr);
    ierr = PetscOptionsSetValue( "-pc_type", "lu" );  CHKERRQ(ierr);
    ierr = PetscOptionsSetValue( "-ksp_monitor", " " );  CHKERRQ(ierr);
    ierr = PetscOptionsSetValue( "-log_info", " " );  CHKERRQ(ierr);
  }
  
  return 0;
}

//PETSc_LinSysCore::setLookup=================================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setLookup"
int PETSc_LinSysCore::setLookup( Lookup& lookup ) 
{
  //lookup_ = &lookup;

  return 0;
}

//PETSc_LinSysCore::setGlobalOffsets======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setGlobalOffsets"
int PETSc_LinSysCore::setGlobalOffsets( int len, int* nodeOffsets,
					int* eqnOffsets, int* blkEqnOffsets ) 
{
  int ii,ierr;
  PRINT3(nodeOffsets[len-1],eqnOffsets[len-1],blkEqnOffsets[len-1]);

  if( proc_gnode_ != NULL ) PetscFree(proc_gnode_);
  ierr = PetscMalloc(sizeof(int)*len,&proc_gnode_); CHKERRQ(ierr);
  for( ii=0 ; ii<len ; ii++ ) proc_gnode_[ii] = nodeOffsets[ii];
  if(numProcs_+1!=len){PRINT2(numProcs_,len-1);ERR_RETURN("np+1!=len");}

  if( proc_globalEq_ != NULL ) PetscFree(proc_globalEq_);
  ierr = PetscMalloc(sizeof(int)*len,&proc_globalEq_);CHKERRQ(ierr);
  for( ii=0 ; ii<len ; ii++ ) proc_globalEq_[ii] = eqnOffsets[ii];

  localStartRow_ = eqnOffsets[thisProc_];
  localEndRow_ = eqnOffsets[thisProc_+1] - 1;
  numLocalRows_ = localEndRow_ - localStartRow_ + 1;
  numGlobalRows_ = eqnOffsets[numProcs_];

  if( proc_lmEq_ != NULL ) PetscFree(proc_lmEq_);
  ierr = PetscMalloc(sizeof(int)*len,&proc_lmEq_);CHKERRQ(ierr);
  for(ii=0 ; ii<len ; ii++) proc_lmEq_[ii] = blkEqnOffsets[ii]-nodeOffsets[ii];

  localStartRowLm_ = proc_lmEq_[thisProc_];
  localEndRowLm_ = proc_lmEq_[thisProc_+1] - 1;
  numLocalRowsLm_ = localEndRowLm_ - localStartRowLm_ + 1;
  numGlobalRowsLm_ = proc_lmEq_[numProcs_];

  if( proc_primEq_ != NULL ) PetscFree(proc_primEq_);
  ierr = PetscMalloc(sizeof(int)*len,&proc_primEq_);CHKERRQ(ierr);
  for( ii=0 ; ii<len ; ii++ ) proc_primEq_[ii] = eqnOffsets[ii]-proc_lmEq_[ii];

  localStartRowPrim_ = proc_primEq_[thisProc_];
  localEndRowPrim_ = proc_primEq_[thisProc_+1] - 1;
  numLocalRowsPrim_ = localEndRowPrim_ - localStartRowPrim_ + 1;
  numGlobalRowsPrim_ = proc_primEq_[numProcs_];

  return 0;
}

//PETSc_LinSysCore::putNodalFieldData========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::putNodalFieldData"
int PETSc_LinSysCore::putNodalFieldData( int fieldID, 
					 int fieldSize,
					 int* nodeNumbers, 
					 int numNodes,
					 const double* data)
{
  PRINT3(fieldID,fieldSize,numNodes);
  return 0; 
}

//PETSc_LinSysCore::setConnectivities=========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setConnectivities"
int PETSc_LinSysCore::setConnectivities(GlobalID elemBlock,
					 int numElements,
					 int numNodesPerElem,
					 const GlobalID* elemIDs,
					 const int* const* connNodes)
{
  return 0; 
}

//PETSc_LinSysCore::setMultCREqns======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::setMultCREqns"
int PETSc_LinSysCore::setMultCREqns(int multCRSetID,
				    int numCRs, 
				    int numNodesPerCR,
				    int** nodeNumbers, 
				    int** eqnNumbers,
				    int* fieldIDs,
				    int* multiplierEqnNumbers )
{
  PRINT3(multCRSetID,numCRs,numNodesPerCR); 
  int ii,jj;
  
  // keep track for regularization!!!!
  for( ii = 0 ; ii < numCRs ; ii++ ) { 
    PRINT1( multiplierEqnNumbers[ii] );
    for( jj = 0 ; jj < numNodesPerCR ; jj++ ){
      PRINT2( nodeNumbers[ii][jj], eqnNumbers[ii][jj] );
      // prometheus -- add into exta_edges list
    }
  }
  for( jj = 0 ; jj < numNodesPerCR ; jj++ ){
    PRINT1( fieldIDs[jj] );
  }

  // keep track of count
  if( numNodesPerCR > maxNodesPerCR_ ) maxNodesPerCR_ = numNodesPerCR;
  
  return 0; 
}

//PETSc_LinSysCore::setPenCREqns========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setPenCREqns"
int PETSc_LinSysCore::setPenCREqns(int penCRSetID,
				   int numCRs, 
				   int numNodesPerCR,
				   int** nodeNumbers, 
				   int** eqnNumbers,
				   int* fieldIDs)
{
  PRINT2( penCRSetID, numCRs );
  return 0;
}

//PETSc_LinSysCore::setMatrixStructure========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setMatrixStructure"
int PETSc_LinSysCore::setMatrixStructure( int** ptColIndices, 
					  int* ptRowLengths,
					  int** blkColIndices, 
					  int* blkRowLengths,
					  int* ptRowsPerBlkRow )
{
  //
  //This function is where we establish the structures/objects associated
  //with the linear algebra library. i.e., do initial allocations, etc.
  //
  int ierr,ii,kk,tt,geq,myeq,*nnz=NULL,*offnnz=NULL,*Cnnz=NULL,*Coffnnz=NULL;
  int lamb,proc;
  Mat A, C = NULL; 
  Vec vecx, vecp = NULL;
  const int nlocalNd = proc_gnode_[thisProc_+1] - proc_gnode_[thisProc_];
  PRINT;

  if(proc_globalEq_ == NULL){std::cerr << "proc_globalEq_ == NULL"<<std::endl;return(-1);}  
  if(proc_globalEq_[numProcs_] != numGlobalRows_){std::cerr<<"ng!=pg[np]"<<std::endl;return(-1);}  
  if( K_ != NULL ) { 
    delete K_;  K_ = NULL; 
    delete x_; x_ = NULL; delete init_guess_; init_guess_ = NULL;
    for(ii=0; ii < numRHSs_; ii++){if(b_!=NULL){delete b_[ii]; b_[ii]=NULL;}}
    PetscFree( b_ ); b_ = NULL;
  }
  // get ndf_ (need to change for prometheus)
  if( ndf_ != 0 ) { std::cerr<<"ndf_ != 0"<<std::endl;return(-1); }  
  
  if( nlocalNd ) { 
    for( ii = 1, ndf_ = ptRowsPerBlkRow[0] ; ii < nlocalNd ; ii++ ) {
      if( ptRowsPerBlkRow[ii] != ndf_ ) { ndf_ = 1; break; }
    }
  }
  else ndf_ = 9999;
  MPI_Allreduce( &ndf_, &tt, 1, MPI_INT, MPI_MIN, FEI_PETSc_Comm_ );
  if( !nlocalNd ) ndf_ = tt;

  if( tt > 1 ) {
    MPI_Allreduce( &ndf_, &kk, 1, MPI_INT, MPI_MAX, FEI_PETSc_Comm_ );
    if( tt != kk ) ndf_ = 1;
    else if( tt != ndf_ || kk != ndf_ ){std::cerr<<"wrong ndf_"<<std::endl;return(-1);}
  }
  else ndf_ = 1;
  PRINT1(ndf_);

  // stiffness matrix  -- need to call prometheus init here and take mat,vec
  if( numLocalRowsLm_ ){
    ierr = PetscMalloc(numLocalRowsLm_*sizeof(int),&Cnnz);CHKERRQ(ierr);
    ierr = PetscMemzero(Cnnz,numLocalRowsLm_*sizeof(int)); CHKERRQ(ierr);
    ierr = PetscMalloc(numLocalRowsLm_*sizeof(int),&Coffnnz);CHKERRQ(ierr);
    ierr = PetscMemzero(Coffnnz,numLocalRowsLm_*sizeof(int)); CHKERRQ(ierr);
  }
  if( ndf_ == 1 ) {
    if( numLocalRowsPrim_ ) {
      ierr = PetscMalloc(numLocalRowsPrim_*sizeof(int),&nnz);CHKERRQ(ierr); 
      ierr = PetscMalloc(numLocalRowsPrim_*sizeof(int),&offnnz);CHKERRQ(ierr);
    }
    for( ii = 0 ; ii < nlocalNd ; ii++ ) {
      const int *pi = ptColIndices[ii], len = ptRowLengths[ii];
      for( kk = nnz[ii] = offnnz[ii] = 0 ; kk < len ; kk++ ) {
	geq = pi[kk]; 
	if(geq >= numGlobalRows_){std::cerr<<"geq >= numGlobalEqns"<<std::endl;return(-1);}
	ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
	if( !lamb ) {
	  if( proc != thisProc_ ) offnnz[ii]++;
	  else nnz[ii]++;
	}
	else {  
	  tt = myeq - localStartRowLm_;
	  if( tt < 0 || tt >= numLocalRowsLm_ ) {
	    // Send 'proc' *row* ('myeq') and *column* ii+localStartRowPrim_
	    // recv and add to Cnnz and add: 
	    // clist = row_colList.Find(myeq+1); clist->AddTail(ii+lp+1);
	  }  
	  else{
	    //Cnnz[tt]++; 
	    // add: clist = row_colList.Find(myeq+1); clist->AddTail(ii+lp+1);
	  }
	} 
      }
    }
    // send/recv lists, add to Cnnz and Coffnnz
    // regulariation to nnz,offnz
    ierr = MatCreateMPIAIJ( FEI_PETSc_Comm_, numLocalRowsPrim_, 
			    numLocalRowsPrim_, numGlobalRowsPrim_, 
			    numGlobalRowsPrim_, 0, nnz, 0, offnnz, &A );
    CHKERRQ(ierr); 
    if( nnz ) { PetscFree( nnz );   PetscFree( offnnz ); }
  }
  else {
    if( nlocalNd ) {
      ierr = PetscMalloc(nlocalNd*sizeof(int),&nnz); CHKERRQ(ierr);
      ierr = PetscMalloc(nlocalNd*sizeof(int),&offnnz);CHKERRQ(ierr);
    }
    for( ii = 0 ; ii < nlocalNd ; ii++ ) {
      const int *pi = blkColIndices[ii], len = blkRowLengths[ii];
      for( kk = nnz[ii] = offnnz[ii] = 0 ; kk < len ; kk += ndf_ ) {
	geq = pi[kk]; 
	if(geq >= numGlobalRows_){std::cerr<<"geq >= numGlobalEqns"<<std::endl;return(-1);}
	ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
	if( !lamb ) {
	  if( proc != thisProc_ ) offnnz[ii]++;
	  else nnz[ii]++;
	} 
	// else  ****from above*****
      }
      offnnz[ii] *= ndf_;
      nnz[ii] *= ndf_;
    } 
    // send/recv lists, add to Cnnz and Coffnnz
    // and regulariation to nnz,offnnz
    ierr = MatCreateMPIBAIJ( FEI_PETSc_Comm_, ndf_, numLocalRowsPrim_, 
			     numLocalRowsPrim_, numGlobalRowsPrim_, 
			     numGlobalRowsPrim_, 0, nnz, 
			     0, offnnz, &A );
    CHKERRQ(ierr);
    if( nnz ) { PetscFree( nnz );   PetscFree( offnnz ); }
  } 
  ierr = MatSetOption( A, MAT_KEEP_ZEROED_ROWS ); CHKERRQ(ierr);  
  //ierr = MatSetOption( A, MAT_NEW_NONZERO_ALLOCATION_ERR );  CHKERRQ(ierr);
  // constraint eqs
  if( numGlobalRowsLm_ == 0 ) C = NULL;
  else { 
    // LM matrix C
    if( maxNodesPerCR_ == 0 ) maxNodesPerCR_ = 1; // fix in new FEI...
    tt = maxNodesPerCR_ * ndf_; // estimate of max nz per row
    ierr = MatCreateMPIAIJ( FEI_PETSc_Comm_, 
			    numLocalRowsLm_, numLocalRowsPrim_, 
			    numGlobalRowsLm_, numGlobalRowsPrim_, 
			    tt, NULL, 0, NULL, &C );
    CHKERRQ(ierr);
    //ierr = MatSetOption( C, MAT_NEW_NONZERO_ALLOCATION_ERR );  CHKERRQ(ierr);
    
    // normilization vector (diagonal matrix)
    Vec PenStiff = NULL; 
    ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsLm_, numGlobalRowsLm_, 
			 &PenStiff );
    CHKERRQ(ierr)
    Pen_stiff_void_ = PenStiff;
    ierr = VecSet( &zero, PenStiff );  CHKERRQ(ierr);  
    if( Cnnz ) { PetscFree( Cnnz );   PetscFree( Coffnnz ); }
  }
  K_ = new PromCRMat( A, C );
  
  // b array
  if( currentRHS_ < 0 ) currentRHS_ = 0;
  if( numRHSs_ < 0 ) { std::cerr<< "numRHSs_ < 0" <<std::endl;return(-1); }
  ierr = PetscMalloc( numRHSs_*sizeof(PromCRVec*),&b_);CHKERRQ(ierr);
  for( ii = 0 ; ii < numRHSs_ ; ii++ ) b_[ii] = NULL;
  // Vecs
  for( ii = 0; ii < numRHSs_; ii++ ){
    ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsPrim_,
			 numGlobalRowsPrim_, &vecx); CHKERRQ(ierr);  
    ierr = VecSet( &zero, vecx );  CHKERRQ(ierr);  
    if( C != NULL ) { // LM Vecs
      ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsLm_, numGlobalRowsLm_,
			   &vecp); CHKERRQ(ierr);  
      ierr = VecSet( &zero, vecp );  CHKERRQ(ierr);
    }
    b_[ii] = new PromCRVec(vecx,vecp); 
  }
  //
  ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsPrim_, numGlobalRowsPrim_,
 		       &vecx); CHKERRQ(ierr);  
  ierr = VecSet( &zero, vecx );  CHKERRQ(ierr);
  if( C != NULL ) { // LM Vecs
    ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsLm_, numGlobalRowsLm_,
			 &vecp); CHKERRQ(ierr); 
    ierr = VecSet( &zero, vecp );  CHKERRQ(ierr);
  }
  init_guess_ = new PromCRVec(vecx,vecp); 
  //
  ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsPrim_, numGlobalRowsPrim_,
		       &vecx); CHKERRQ(ierr); 
  ierr = VecSet( &zero, vecx );  CHKERRQ(ierr);
  if( C != NULL ) { // LM Vecs
    ierr = VecCreateMPI( FEI_PETSc_Comm_, numLocalRowsLm_, numGlobalRowsLm_,
			 &vecp ); CHKERRQ(ierr); 
    ierr = VecSet( &zero, vecp );  CHKERRQ(ierr);
  }
  x_ = new PromCRVec(vecx,vecp);  
  
  return 0;
}

//PETSc_LinSysCore::resetMatrixAndVector======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::resetMatrixAndVector"
int PETSc_LinSysCore::resetMatrixAndVector( double s ) 
{
  PRINTD(s); 
  int ierr;

  ierr = resetMatrix(s); CHKERRQ(ierr); 
  ierr = resetRHSVector(s); CHKERRQ(ierr); 

  return 0;
}

//PETSc_LinSysCore::resetMatrix======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::resetMatrix"
int PETSc_LinSysCore::resetMatrix( double s ) 
{
  PRINTD(s); 
  if( s != 0.0 ){ std::cerr<< "s != 0." <<std::endl;return(-1); }
  
  if( K_ != NULL ) { int ierr = K_->ZeroEntries(); CHKERRQ(ierr); }

  dirty_system_ = TRUE;
  return 0;
}

//PETSc_LinSysCore::resetRHSVector======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::resetRHSVector"
int PETSc_LinSysCore::resetRHSVector( double s ) 
{
  PRINTD(s); 
  //
  if( b_ != NULL ) {
    PromCRVec *bb = b_[currentRHS_];
    if( bb == NULL ) { std::cerr<<"bb == NULL"<<std::endl;return(-1); }
    int ierr = bb->Set( s ); CHKERRQ(ierr);
  }
  return 0;
}

//PETSc_LinSysCore::putIntoSystemMatrix=======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::putIntoSystemMatrix(point)"
int PETSc_LinSysCore::putIntoSystemMatrix(int numPtRows, 
					  const int* ptRows,
					  int numPtCols, 
					  const int* ptCols,
					  const double* const* values)
{
  
  return sumIntoSystemMatrix_private( numPtRows, ptRows, numPtCols, ptCols,
				      values, (int)INSERT_VALUES );
}

//PETSc_LinSysCore::sumIntoSystemMatrix=======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::sumIntoSystemMatrix(point)"
int PETSc_LinSysCore::sumIntoSystemMatrix(int numPtRows, 
					  const int* ptRows,
					  int numPtCols, 
					  const int* ptCols,
					  const double* const* values)
{
  return sumIntoSystemMatrix_private( numPtRows, ptRows, numPtCols, ptCols,
				      values, (int)ADD_VALUES );
}

//PETSc_LinSysCore::sumIntoSystemMatrix_private=======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::sumIntoSystemMatrix_private(point)"
int PETSc_LinSysCore::sumIntoSystemMatrix_private(int numPtRows, 
						  const int* ptRows,
						  int numPtCols, 
						  const int* ptCols,
						  const double* const* values, int add_type )
{
  //PRINT2(numPtRows,numPtCols);
  if ( K_ == NULL ){std::cerr<<"K_ == NULL"<<std::endl;return(-1);}
  int ii,jj,ierr,geq,myeq,lamb,proc,nc,lambI=99;
  Mat A = K_->A_, C = K_->C_;
  
  if (numPtRows <= 0 || numPtCols <= 0) {
    FEI_CERR << "PETSc_LinSysCore::sumIntoRow: numPtCols: " << numPtCols << FEI_ENDL;
    FEI_CERR << "PETSc_LinSysCore::sumIntoRow: numPtRows: " << numPtRows << FEI_ENDL;
    return(0);
  }
  int *idsI,*idsJ,*maskJ;
  double *vv,*pv;
  ierr = PetscMalloc(numPtRows*sizeof(int),&idsI);CHKERRQ(ierr);
  ierr = PetscMalloc(numPtCols*sizeof(int),&idsJ);CHKERRQ(ierr);
  ierr = PetscMalloc(numPtCols*sizeof(int),&maskJ);CHKERRQ(ierr);
  ierr = PetscMalloc(numPtRows*numPtCols*sizeof(double),&vv);CHKERRQ(ierr);
  pv = vv;

  // cols
  for( jj = 0, nc = 0 ; jj < numPtCols ; jj++ ){
    geq = ptCols[jj];
    ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
    if( lamb ) maskJ[jj] = 0; // skip C^T
    else {
      idsJ[nc++] = myeq;
      maskJ[jj] = 1;
    }
  }
  if( nc ) { // else pure C^T 
    // rows
    for( ii = 0 ; ii < numPtRows ; ii++ ){
      geq = ptRows[ii];
      if( geq > localEndRow_ ){ std::cerr<<"geq > localEndRow_"<<std::endl;return(-1);}
      if( geq < localStartRow_ ){ std::cerr<<"geq < localStartRow_"<<std::endl;return(-1);}
      ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
      if( lambI == 99 ) lambI = lamb; // get lamb
      else if( lamb != lambI ){std::cerr<<"mixed primary and dual rows"<<std::endl;return(-1);}
      idsI[ii] = myeq;
      const double *varr = values[ii];
      for( jj = 0 ; jj < numPtCols ; jj++ ) if( maskJ[jj] ) *pv++ = varr[jj];
    }
    if(pv-vv!=numPtRows*nc){std::cerr<<"numPtRows*nc"<<std::endl;return(-1);}
    // add v
    if( lambI ){
      ierr = MatSetValues(C,numPtRows,idsI,nc,idsJ,vv,INSERT_VALUES);
      CHKERRQ(ierr);
    }
    else {
      ierr = MatSetValues(A,numPtRows,idsI,nc,idsJ,vv,(InsertMode)add_type);
      CHKERRQ(ierr);
      dirty_system_ = TRUE;
    }
  }
  PetscFree(idsI); PetscFree(idsJ); PetscFree(vv); PetscFree(maskJ);
  
  return 0;
}

//PETSc_LinSysCore::getMatrixRowLength=======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::getMatrixRowLength"
int PETSc_LinSysCore::getMatrixRowLength(int row, int& length)
{
  Mat A = K_->A_, C = K_->C_; 
  double *matrowv;
  int ierr,ncols,*colids,myeqi,proc,lamb;
  
  ierr = eq_map_private( row, myeqi, lamb, proc ); CHKERRQ(ierr);
  if( lamb ) {
    ierr = MatGetRow( C, myeqi, &ncols, &colids, &matrowv );  CHKERRQ(ierr);
  }
  else {
    ierr = MatGetRow( A, myeqi, &ncols, &colids, &matrowv );  CHKERRQ(ierr);
  }

  // ret
  length = ncols; 
    
  if( lamb ) {
    ierr = MatRestoreRow(C,myeqi, &ncols, &colids, &matrowv); CHKERRQ(ierr);
  }
  else {
    ierr = MatRestoreRow(A,myeqi, &ncols, &colids, &matrowv); CHKERRQ(ierr);
  }

  return(0);
}

//PETSc_LinSysCore::getMatrixRow=======================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::getMatrixRow"
int PETSc_LinSysCore::getMatrixRow(int row, double* coefs,
				   int* indices,
				   int len, int& rowLength)
{
  Mat A = K_->A_, C = K_->C_; 
  double *matrowv;
  int ierr,ii,ncols,*colids,myeqi,proc,lamb;
  
  ierr = eq_map_private( row, myeqi, lamb, proc ); CHKERRQ(ierr);
  if( lamb ) {
    ierr = MatGetRow( C, myeqi, &ncols, &colids, &matrowv );  CHKERRQ(ierr);
  }
  else {
    ierr = MatGetRow( A, myeqi, &ncols, &colids, &matrowv );  CHKERRQ(ierr);
  }

  static char seterrq2_string[] = "provide ros len(%d) < %d";

  // ret
  rowLength = ncols; 
  if(len<ncols)SETERRQ2(1,seterrq2_string,len,ncols);
  for(ii=0;ii<ncols;ii++){indices[ii] = colids[ii]; coefs[ii] = matrowv[ii];}
  
  if( lamb ) {
    ierr = MatRestoreRow(C,myeqi, &ncols, &colids, &matrowv); CHKERRQ(ierr);
  }
  else {
    ierr = MatRestoreRow(A,myeqi, &ncols, &colids, &matrowv); CHKERRQ(ierr);
  }

  return(0);
}

//PETSc_LinSysCore::sumIntoSystemMatrix==================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::sumIntoSystemMatrix(blocked)"
int PETSc_LinSysCore::sumIntoSystemMatrix(int numPtRows,
					  const int* ptRows,
					  int numPtCols,
					  const int* ptCols,
					  int numBlkRows,
					  const int* blkRows,
					  int numBlkCols,
					  const int* blkCols,
					  const double* const* values ) 
{
  if( K_ == NULL ){ std::cerr<< "K_ == NULL" <<std::endl;return(-1); }
  int ii,jj,kk,geq,myeq,lamb,lambI=99,ierr,proc,nbc;

  if( ndf_ == 1 ) { 
    ierr = sumIntoSystemMatrix( numPtRows, ptRows, numPtCols, ptCols, values );
    CHKERRQ(ierr);
  }
  else {
    const int nv = numPtRows*numPtCols, matndf = numPtRows/numBlkRows;
    Mat A = K_->A_, C = K_->C_;
    double *pv,*vv;
    int *rows,*cols,*maskJ;
    ierr = PetscMalloc(numBlkRows*sizeof(int),&rows);CHKERRQ(ierr);
    ierr = PetscMalloc(numBlkCols*sizeof(int),&cols);CHKERRQ(ierr);
    ierr = PetscMalloc(numBlkCols*sizeof(int),&maskJ);CHKERRQ(ierr);
    ierr = PetscMalloc(nv*sizeof(double),&vv);CHKERRQ(ierr);
    ierr = MatGetBlockSize( A, &ii );  CHKERRQ(ierr);
    if( matndf != ii ){ std::cerr<<"block size messes up???"<<std::endl;return(-1);}

    // cols
    for( jj = 0, nbc = 0, kk = 0 ; kk < numBlkCols ; jj += matndf, kk++ ){
      geq = ptCols[jj];
      ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
      if( lamb ) maskJ[kk] = 0; // skip C^T
      else {
	cols[nbc++] = myeq/matndf;
	maskJ[kk] = 1;
      }
    }
    if( kk != numBlkCols ){std::cerr<<"kk != numBlkCols"<<std::endl;return(-1);}
    if( nbc > numBlkCols ){ std::cerr<<"nbc > numBlkCols"<<std::endl;return(-1); }
    if( nbc ) { // else pure C^T 
      // rows
      for( ii = 0, kk = 0, pv = vv ; kk < numBlkRows ; ii += matndf, kk++ ){
	geq = ptRows[ii];
	if( geq > localEndRow_ ){ERR_RETURN("geq > localEndRow_");}
	if( geq < localStartRow_ ){ERR_RETURN("geq < localStartRow_");}
	ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
	if( lambI == 99 ) lambI = lamb; // get lamb
	else if( lamb != lambI ){ERR_RETURN("mixed primary and dual rows");}
	rows[kk] = myeq/matndf;
	for( int iii = 0, jjj ; iii < matndf ; iii++ ) {
	  const double *varr = values[ii + iii];
	  for( jj = jjj = 0 ; jj < numBlkCols ; jj++, jjj += matndf ) {
	    if( maskJ[jj] ) for( int j=0; j<matndf ; j++ ) *pv++ = varr[jjj+j];
	  }
	}
      }
      if( kk != numBlkRows ){ERR_RETURN("kk != numBlkRows");}
      if(pv-vv!=numPtRows*nbc*matndf){ERR_RETURN("numPtRows*nbc");}
      // add
      if( lambI ) {
	ierr = MatSetValues(C,numPtRows,rows,numPtCols,cols,vv,INSERT_VALUES);
	CHKERRQ(ierr);
      }
      else{
	ierr = MatSetValuesBlocked(A,numBlkRows,rows,nbc,cols,vv,ADD_VALUES);
	CHKERRQ(ierr);
	dirty_system_ = TRUE;
      }
    }
    PetscFree(vv); PetscFree(rows); PetscFree(cols);  PetscFree(maskJ); 
  }

  return 0;
}

//PETSc_LinSysCore::putInitialGuess===========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::putInitialGuess"
int PETSc_LinSysCore::putInitialGuess(const int* eqnNumbers,
				       const double* values, int len) 
{
  //
  //This function scatters (puts) values into the linear-system's soln vector.
  //
  // num is how many values are being passed,
  // indices holds the global 'row-numbers' into which the values go,
  // and values holds the actual coefficients to be scattered.
  //
  PRINT1(len);

  if( len ) {
    int ii,jj,geq,myeq,lamb,lamb1=99,ierr,proc,*rows; double *vv;
    Vec vecx = init_guess_->x_, vecp = init_guess_->p_;
    ierr = PetscMalloc(len*sizeof(int),&rows);CHKERRQ(ierr);
    ierr = PetscMalloc(len*sizeof(double),&vv);CHKERRQ(ierr);
    for( ii = jj = 0 ; ii < len ; ii++ ){
      geq = eqnNumbers[ii];
      ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
      if( proc == thisProc_){ vv[jj] = values[ii]; rows[jj] = myeq; jj++; }
      if( lamb1 == 99 ) lamb1 = lamb;
      else if( lamb1 != lamb ) { ERR_RETURN( "lamb1 != lamb" ); }
    }
    if( jj ) {
      if( lamb1 ) {
	ierr = VecSetValues( vecp, jj, rows, vv, INSERT_VALUES); CHKERRQ(ierr);
      }
      else {
	ierr = VecSetValues( vecx, jj, rows, vv, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
    PetscFree(rows); PetscFree(vv); 
    dirty_guess_ = TRUE;
  }

  return 0;
}
 
//PETSc_LinSysCore::getFromRHSVector==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getFromRHSVector"
int PETSc_LinSysCore::getFromRHSVector(int num, double* values,
				       const int* indices) 
{
  if(num<=0) return 0;
  int ii,geq,myeq,lamb1=99,proc,ierr,lamb,*inds,inds_buff[1024];
  PromCRVec *b_cr = b_[currentRHS_];
  if(b_cr == NULL){ ERR_RETURN( "b_[currentRHS_] == NULL" ); }  
  
  if( num > 1024 ){ 
    ierr = PetscMalloc(num*sizeof(int),&inds); CHKERRQ(ierr);
  }
  else inds = inds_buff;
  
  for( ii = 0 ; ii < num ; ii++ ){ 
    geq = indices[ii]; 
    ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
    static char seterrq_str1[] = "pm1!=pm";
    if( lamb1 == 99 ) lamb1 = lamb; 
    else if(lamb1!=lamb){SETERRQ(1,seterrq_str1);}
    inds[ii] = myeq; // overwrites 'indices' with ly indices
  }
  if( lamb1 ){ 
    const int lm0 = proc_lmEq_[thisProc_]; double *pv;
    ierr = VecGetArray( b_cr->p_, &pv );    CHKERRQ(ierr);
    for(int ii=0;ii<num;ii++) values[ii] = pv[inds[ii] - lm0];
    ierr = VecRestoreArray( b_cr->p_, &pv );    CHKERRQ(ierr);
  }
  else{
    const int prim0 = proc_primEq_[thisProc_]; double *pv;
    ierr = VecGetArray( b_cr->x_, &pv );    CHKERRQ(ierr);
    for(int ii=0;ii<num;ii++) values[ii] = pv[inds[ii] - prim0];
    ierr = VecRestoreArray( b_cr->x_, &pv );    CHKERRQ(ierr);
  } 
  if( inds != inds_buff ) PetscFree(inds); 
  return 0;
}

//PETSc_LinSysCore::putIntoRHSVector==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::putIntoRHSVector"
int PETSc_LinSysCore::putIntoRHSVector(int num, const double* values,
				       const int* indices) 
{
  return sumIntoRHSVector_private( num, values, indices, (int)INSERT_VALUES) ;
}

//PETSc_LinSysCore::sumIntoRHSVector==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::sumIntoRHSVector"
int PETSc_LinSysCore::sumIntoRHSVector(int num, const double* values,
				       const int* indices) 
{
  return sumIntoRHSVector_private( num, values, indices, (int)ADD_VALUES) ;
}

//PETSc_LinSysCore::sumIntoRHSVector_private==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::sumIntoRHSVector_private"
int PETSc_LinSysCore::sumIntoRHSVector_private(int num, const double* values,
					       const int* indices, int add_type) 
{
  //
  //This function scatters (accumulates) values into the linear-system's
  //currently selected RHS vector.
  //
  // num is how many values are being passed,
  // indices holds the global 'row-numbers' into which the values go,
  // and values holds the actual coefficients to be scattered.
  // 
  //PRINT1(num);  
  PromCRVec *bb = b_[currentRHS_]; Vec vec;
  int ierr,ii,geq,myeq,lamb,lamb1=99,proc,*rows,__buff[1024];
  if(bb == NULL){ ERR_RETURN( "b_[currentRHS_] == NULL" ); }
  
  if(num> 1024){ierr = PetscMalloc(num*sizeof(int),&rows);CHKERRQ(ierr);}
  else rows = __buff;
  for(ii=0;ii<num;ii++){ 
    geq = indices[ii]; 
    ierr = eq_map_private( geq, myeq, lamb, proc ); CHKERRQ(ierr);
    if( lamb1 == 99 ) lamb1 = lamb; 
    else if(lamb1!=lamb||proc!=thisProc_){ERR_RETURN("pm1!=pm || pe!=mype");}
    rows[ii] = myeq;
  }
  
  vec = lamb1 ? bb->p_ : bb->x_;
  if( INSERT_VALUES == add_type ) {
    double *vv,*vv2,__fbuff[1024];
    if(num> 1024){ierr = PetscMalloc(num*sizeof(double),&vv);CHKERRQ(ierr);}
    else vv = __fbuff;
    const int my0 = lamb1 ? proc_lmEq_[thisProc_] : proc_primEq_[thisProc_];
    //const int myNext0 = lamb1 ? proc_lmEq_[mype_+1]:proc_primEq_[mype_+1];
    ierr = VecGetArray( vec, &vv2 );    CHKERRQ(ierr);
    for( int ii = 0 ; ii < num ; ii++ ){
      //assert( rows[ii] >=  my0 && rows[ii] < myNext0 );
      vv[ii] = values[ii] - vv2[ rows[ii] - my0 ];
    }
    ierr = VecRestoreArray( vec, &vv2 ); CHKERRQ(ierr); 
    ierr = VecSetValues( vec, num, rows, vv, ADD_VALUES ); CHKERRQ(ierr);
    if( vv != __fbuff) PetscFree(vv); 
  }
  else {
    ierr = VecSetValues( vec, num, rows, values, ADD_VALUES ); CHKERRQ(ierr);
  } 

  if( rows != __buff) PetscFree(rows); 
  return 0;
}

//PETSc_LinSysCore::matrixLoadComplete========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::matrixLoadComplete"
int PETSc_LinSysCore::matrixLoadComplete() 
{
  PRINT;
  int ierr; Mat A = K_->A_, C = K_->C_;
  PromCRVec *bb = b_[currentRHS_];
  // 
  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
  if( C != NULL ) {
    ierr = MatAssemblyBegin( C, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
    ierr = MatAssemblyEnd( C, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
    ierr = VecAssemblyBegin( bb->p_ );                     CHKERRQ(ierr);
    ierr = VecAssemblyEnd( bb->p_ );      CHKERRQ(ierr);
  }
  //
  ierr = VecAssemblyBegin( bb->x_ );                     CHKERRQ(ierr);
  ierr = VecAssemblyEnd( bb->x_ );      CHKERRQ(ierr);
  // call prometheus build grids
  return 0;
}

//PETSc_LinSysCore::enforceEssentialBC========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::enforceEssentialBC"
int PETSc_LinSysCore::enforceEssentialBC(int* globalEqn,
					  double* alpha,
					  double* gamma, int len) 
{
  if( len ) {
    PETSc_EssBCData *bc = new PETSc_EssBCData(globalEqn,alpha,gamma,len);
    bc->next = essbcdata_;
    essbcdata_ = bc;
  }
  return 0;
}

//PETSc_LinSysCore::enforceEssentialBC_private ===============================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::enforceEssentialBC_private"
int PETSc_LinSysCore::enforceEssentialBC_private( int* globalEqn, 
						  double* alpha, // in/out!!!
						  double* gamma, int len ) 
{
  //
  //This function must enforce an essential boundary condition on each local
  //equation in 'globalEqn'. This means, that the following modification
  //should be made to A and b, for each globalEqn:
  //
  //for(each local equation i){
  //   for(each column j in row i) {
  //      if (i==j) b[i] = gamma/alpha;
  //      else b[j] -= (gamma/alpha)*A[j,i];
  //   }
  //}
  //
  //all of row 'globalEqn' and column 'globalEqn' in A should be zeroed,
  //except for 1.0 on the diagonal.
  //
  PRINT1(len); 

  double *matrowv, *bvec; 
  int ierr,self,jj,ii,kk,*colids,*bids,ncols,nloc=0,proc; 
  int *localrows = NULL,*sizes=NULL,**idsarr=NULL;
  PromCRVec *bb = b_[currentRHS_]; Vec b = bb->x_;
  Mat A = K_->A_;
  
  ierr = PetscMalloc(len*sizeof(int),&localrows);CHKERRQ(ierr);
  ierr = PetscMalloc(len*sizeof(int),&sizes);CHKERRQ(ierr);
  ierr = PetscMalloc(len*sizeof(int*),&idsarr);CHKERRQ(ierr);
  for( ii=0; ii<len; ii++) {
    //if globalEqn[i] is local, we'll diagonalize the row and column.
    int colj, myeqi, lamb, geqi = globalEqn[ii];
    if( geqi >= localStartRow_ && geqi <= localEndRow_ ){
      const double rhs_term = gamma[ii]/alpha[ii];
      alpha[ii] = rhs_term; // output!!!
      ierr = eq_map_private( geqi, myeqi, lamb, proc ); CHKERRQ(ierr);
      if( lamb || proc!=thisProc_){ ERR_RETURN("lamb || pe!=mype");}
      // symmetric matrix only!!!, should use PETSc MatGetColumnVector
      ierr = MatGetRow( A, myeqi, &ncols, &colids, &matrowv );
      CHKERRQ(ierr);
      ierr = PetscMalloc(ncols*sizeof(double),&bvec);CHKERRQ(ierr);
      ierr = PetscMalloc(ncols*sizeof(int),&bids);CHKERRQ(ierr);
      for( jj = self = kk = 0; jj < ncols; jj++ ) {
	colj = colids[jj];
	if( colj == myeqi ) self++; //diag[nloc] = matrowv[jj];
	else {
	  bids[kk] = colj;
	  bvec[kk++] = -matrowv[jj] * rhs_term; 
	}
      } // end for(jj<rowLength) loop
      // put row back
      ierr = MatRestoreRow(A,myeqi, &ncols, &colids, &matrowv); CHKERRQ(ierr);
      if( kk != ncols-1 ){ ERR_RETURN("ncols-1"); }
      if( self != 1 ){ ERR_RETURN("self!=1"); }
      
      localrows[nloc] = myeqi; // keep for row
      sizes[nloc] = kk;
      idsarr[nloc] = bids;
      nloc++; 
      
      // modify rhs
      ierr = VecSetValues( b, kk, bids, bvec, ADD_VALUES );     CHKERRQ(ierr);
      PetscFree(bvec);       // clean up
    }
  }
  
  // zero rows and cols
  if( nloc ) {  
    PETSc_ZeroEquations *zeqs = new PETSc_ZeroEquations( localrows, sizes, 
							 idsarr, nloc );
    zeqs->next = zeroEqs_;
    zeroEqs_ = zeqs;
  }
  else { 
    PetscFree(localrows); PetscFree(sizes);PetscFree(idsarr);//PetscFree(diag);
  }

  return 0;
}

//PETSc_LinSysCore::enforceOtherBC============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::enforceOtherBC"
int PETSc_LinSysCore::enforceOtherBC( int* globalEqn, double* alpha,
				      double* beta, double* gamma, int len ) 
{
  if( len ) {
    PETSc_OthBCData *bc = new PETSc_OthBCData(globalEqn,alpha,beta,gamma,len);
    bc->next = othbcdata_;
    othbcdata_ = bc;
  }
  return 0;
}

//PETSc_LinSysCore::enforceOtherBC_private ===================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::enforceOtherBC_private"
int PETSc_LinSysCore::enforceOtherBC_private(int* globalEqn, double* alpha,
					      double* beta, double* gamma, 
					      int len) 
{
  //
  //This function must enforce a natural or mixed boundary condition on the
  //equations in 'globalEqn'. This means that the following modification 
  //should be made to A and b:
  //
  //A[globalEqn,globalEqn] += alpha/beta;
  //b[globalEqn] += gamma/beta;
  //
  PRINT1(len); 

  double temp, *bvec = NULL;
  int ierr, myeqi, lamb, geqi,ii,nloc=0,*localrows=NULL,proc;
  PromCRVec *bb = b_[currentRHS_]; Vec b = bb->x_;
  Mat A = K_->A_;

  ierr = PetscMalloc(len*sizeof(int),&localrows);CHKERRQ(ierr);
  ierr = PetscMalloc(len*sizeof(double),&bvec);CHKERRQ(ierr);
  for( ii=0; ii<len; ii++) {
    //if globalEqn[i] is local, we'll diagonalize the row and column.
    geqi = globalEqn[ii];  
    if( geqi >= localStartRow_ && geqi <= localEndRow_){
      ierr = eq_map_private( geqi, myeqi, lamb, proc ); CHKERRQ(ierr);
      if( lamb || proc!=thisProc_ ){ ERR_RETURN("lamb || pe!=mype");}
      // vec
      localrows[nloc] = myeqi; // keep for rows
      bvec[nloc] = gamma[ii]/beta[ii];   // b_[i] += gamma[i]/beta[i];*/
      nloc++;
      // 
      temp = alpha[ii]/beta[ii]; // coefs[j] += alpha[i]/beta[i];
      ierr = MatSetValues( A, 1, &myeqi, 1, &myeqi, &temp, ADD_VALUES );
      CHKERRQ(ierr);
    }
  }
  // fix 
  if( nloc ) {
    ierr = VecSetValues( b, nloc, localrows, bvec, ADD_VALUES );
    CHKERRQ(ierr);
  }
  PetscFree(localrows); PetscFree(bvec);

  return 0;
}

//PETSc_LinSysCore::enforceRemoteEssBCs=======================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::enforceRemoteEssBCs"
int PETSc_LinSysCore::enforceRemoteEssBCs(int numEqns, int* globalEqns,
                                          int** colIndices, int* colIndLen,
                                          double** coefs) 
{
  //
  //globalEqns should hold eqns that are owned locally, but which contain
  //column indices (the ones in colIndices) which are from remote equations
  //on which essential boundary-conditions need to be enforced.
  //
  //This function will only make the modification if the above conditions
  //hold -- i.e., the equation is a locally-owned equation, and the column
  //index is NOT a locally owned equation.
  //
  PRINT; 

  // set prometheus->findgrid0->nodes[gidI-my0].setFixed(dof,ndof_)
  // ** this is done in enforceEssentialBC **
  return 0;
}

//PETSc_LinSysCore::applyBCs_private =========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::applyBCs_private"
int PETSc_LinSysCore::applyBCs_private()
{
  int ierr,ii,jj,mxcol,nloc,nrows;
  PETSc_ZeroEquations *zeqs, *zeqs2;
  PETSc_EssBCData *esbc, *esbc2;
  PETSc_OthBCData *obc, *obc2;
  Mat A = K_->A_, C = K_->C_; 
  int *rows = NULL;double *diags = NULL; IS is = NULL;
  PRINT;

  // OtherBC == calls MatSetValues (diagonal only) -- ADD_VALUES
  obc = othbcdata_ ;
  while( obc ) {
    ierr = enforceOtherBC_private( obc->globalEqn_, obc->alpha_, obc->beta_,
				   obc->gamma_, obc->len_); 
    CHKERRQ(ierr);
    obc2 = obc->next;      delete obc;        obc = obc2; 
  }
  othbcdata_ = NULL;

  // assemble -- collective
  matrixLoadComplete(); 

  // count local rows
  for( nrows = 0, esbc = essbcdata_ ; esbc ; esbc = esbc->next ) {
    for( jj = 0 ; jj < esbc->len_ ; jj++ ) {
      int geqi = esbc->globalEqn_[jj]; 
      if( geqi >= localStartRow_ && geqi <= localEndRow_ )  nrows++;
    }
  }
  const int numrhs = nrows;
  if( numrhs ) {
    ierr = PetscMalloc(numrhs*sizeof(int),&rows);CHKERRQ(ierr);
    ierr = PetscMalloc(numrhs*sizeof(double),&diags);CHKERRQ(ierr);
  }
  esbc = essbcdata_; ii = 0;
  while( esbc ) {
    // calls matgetrow -- makes zeroEqs_
    ierr = enforceEssentialBC_private( esbc->globalEqn_, esbc->alpha_,
				       esbc->gamma_, esbc->len_ );
    CHKERRQ(ierr);
    // collect rhs sets
    for( jj = 0 ; jj < esbc->len_ ; jj++ ) {
      int geqi = esbc->globalEqn_[jj], proc, lamb, myeqi;
      if( geqi >= localStartRow_ && geqi <= localEndRow_ ) {
	ierr = eq_map_private( geqi, myeqi, lamb, proc ); CHKERRQ(ierr);
	if( lamb || proc!=thisProc_ ){ ERR_RETURN("lamb || pe!=mype");}
	rows[ii] = myeqi; diags[ii] = esbc->alpha_[jj]; 
	ii++;
      }
    }
    esbc2 = esbc->next;     delete esbc;        esbc = esbc2;
  }
  essbcdata_ = NULL;
  // set rhs
  PromCRVec *bb = b_[currentRHS_]; Vec b = bb->x_; 
  ierr = VecAssemblyBegin( b );    CHKERRQ(ierr);
  ierr = VecAssemblyEnd( b );      CHKERRQ(ierr);
  if( numrhs ){
    ierr = VecSetValues( b, numrhs, rows, diags, INSERT_VALUES); CHKERRQ(ierr);
    PetscFree(diags);  PetscFree(rows);  diags = NULL; rows = NULL;
  }

  // zero rows and cols
  for( nrows = mxcol = 0, zeqs = zeroEqs_ ; zeqs ; zeqs = zeqs->next ) {
    // get max
    nloc = zeqs->len_;
    for( jj = 0; jj < nloc; jj++ ){
      if( zeqs->sizes_[jj] > mxcol ) mxcol = zeqs->sizes_[jj];
    }
    nrows += nloc;
  }
  // zero rows -- set diag == 1.0  
  if( nrows ) {
    ierr = PetscMalloc(nrows*sizeof(int),&rows);CHKERRQ(ierr);
    ierr = PetscMalloc(nrows*sizeof(double),&diags);CHKERRQ(ierr);
  }  
  for( ii = 0, zeqs = zeroEqs_ ; zeqs ; zeqs = zeqs->next ) {
    nloc = zeqs->len_;
    for( jj = 0; jj < nloc; jj++, ii++ ){
      rows[ii] = zeqs->localrows_[jj]; diags[ii] = 1.0; 
    }
  }
  if( ii != nrows ){ ERR_RETURN("ii != nrows"); }
  // zero rows with old diagonal -- collective
  ierr = ISCreateGeneral( MPI_COMM_SELF, nrows, rows, &is ); CHKERRQ(ierr);
  ierr = MatZeroRows( A, is, diags ); CHKERRQ(ierr);
  if( C != NULL ) {
    Mat transC, tempC; 
    ierr = MatTranspose( C, &transC );      CHKERRQ(ierr);  // transC = C'
    ierr = MatSetOption( transC, MAT_KEEP_ZEROED_ROWS  ); CHKERRQ(ierr);
    ierr = MatZeroRows( transC, is, PETSC_NULL ); CHKERRQ(ierr); // apply BCs
    ierr = MatTranspose( transC, &tempC ); CHKERRQ(ierr); // tempC = C'' (=C)
    ierr = MatCopy( tempC, C, SAME_NONZERO_PATTERN );  CHKERRQ(ierr); // C = tempC
    MatDestroy( transC ) ; MatDestroy( tempC );
  }
  ierr = ISDestroy( is ); CHKERRQ(ierr); 
  if( nrows ){
    PetscFree(diags);  PetscFree(rows);  diags = NULL; rows = NULL;
  }
  // zero cols -- no diagonal
  if( nrows ){ 
    if( zeroEqs_ == NULL ){ ERR_RETURN("zeroEqs_ == NULL"); }
    double *zeros;
    ierr = PetscMalloc(mxcol*sizeof(double),&zeros);CHKERRQ(ierr);
    for( jj = 0; jj < mxcol; jj++ ) zeros[jj] = 0.0;
    zeqs = zeroEqs_;
    while( zeqs ) {
      nloc = zeqs->len_;
      for( jj = 0; jj < nloc; jj++ ) {
	int *ids = zeqs->idsarr_[jj];	
	int sz = zeqs->sizes_[jj];   if(sz > mxcol){ERR_RETURN("sz>mxcol");}
	int col = zeqs->localrows_[jj];
	ierr = MatSetValues( A, sz, ids, 1, &col, zeros, INSERT_VALUES );
	CHKERRQ(ierr);
	PetscFree(ids);
      }
      PetscFree(zeqs->localrows_); PetscFree(zeqs->sizes_); 
      PetscFree(zeqs->idsarr_);
      zeqs2 = zeqs->next;  
      delete zeqs;    
      zeqs = zeqs2;
    }
    PetscFree(zeros); 
    zeroEqs_ = NULL;
  }
  else if (zeroEqs_ != NULL){ERR_RETURN("zeroEqs_ != NULL");}

  // assemble -- collective
  matrixLoadComplete(); 

  return 0;
}

//PETSc_LinSysCore::getMatrixPtr==============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getMatrixPtr"
int PETSc_LinSysCore::getMatrixPtr(Data& data) 
{
  //PRINT1(K_->A_);
  data.setTypeName(matrix_name);
  data.setDataPtr( (void*)K_ );

  return 0;
}


//PETSc_LinSysCore::destroyMatrixData=========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::destroyMatrixData"
int PETSc_LinSysCore::destroyMatrixData( Data& data ) 
{
  PRINT;

  if( strcmp(matrix_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Matrix'.");
  }
  PromCRMat *inA = (PromCRMat*)data.getDataPtr(); 
//  Mat A = inA->A_, C = inA->C_;
  if( inA == K_ ) {
    ERR_RETURN("Dstroying system matrix");
  }

  delete inA;
  data.setDataPtr( NULL );   

  return 0;
}

//PETSc_LinSysCore::copyInMatrix==============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::copyInMatrix"
int PETSc_LinSysCore::copyInMatrix( double scalar, const Data& data ) 
{
  //
  //Overwrites the current internal matrix with a scaled copy of the
  //input argument.
  //
  PRINTD(scalar);
  int ierr;
  if( strcmp(matrix_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Matrix'.");
  }
  PromCRMat *inA = (PromCRMat*)data.getDataPtr(); 
  const Mat A = inA->A_, C = inA->C_; 

  // better have same non-zero pattern!!! -- DIFFERENT_NONZERO_PATTERN
  if( scalar != 1.0 ) {
    Mat    Temp;
    ierr = MatDuplicate( A, MAT_COPY_VALUES, &Temp );  CHKERRQ(ierr);
    ierr = MatScale( &scalar, Temp );  CHKERRQ(ierr);
    ierr = MatCopy( Temp, K_->A_, SAME_NONZERO_PATTERN );  CHKERRQ(ierr);
    ierr = MatDestroy( Temp );  CHKERRQ(ierr);
  }
  else {
    ierr = MatCopy( A, K_->A_, SAME_NONZERO_PATTERN );  CHKERRQ(ierr);
  }

  dirty_system_ = TRUE;

  if( C != NULL ){
    if( scalar != 1.0 ) {
      Mat    Temp;
      ierr = MatDuplicate( C, MAT_COPY_VALUES, &Temp );  CHKERRQ(ierr);
      ierr = MatScale( &scalar, Temp );  CHKERRQ(ierr);
      ierr = MatCopy( Temp, K_->C_, SAME_NONZERO_PATTERN );  CHKERRQ(ierr);
      ierr = MatDestroy( Temp );  CHKERRQ(ierr);
    }
    else {
      ierr = MatCopy( C, K_->C_, SAME_NONZERO_PATTERN );  CHKERRQ(ierr);
    }
  }

  return 0;
}

//PETSc_LinSysCore::copyOutMatrix=============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::copyOutMatrix"
int PETSc_LinSysCore::copyOutMatrix( double scalar, Data& data ) 
{
  //
  //Passes out a scaled (created) copy of the current internal matrix.
  //
  PRINTD(scalar);  
  int ierr; Mat A = K_->A_, C = K_->C_, newA, newC;

  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
  ierr = MatDuplicate( A, MAT_COPY_VALUES, &newA );  CHKERRQ(ierr);
  if( scalar != 1.0 ) {
    ierr = MatScale( &scalar, newA );  CHKERRQ(ierr);
  }

  if( C != NULL ) {
    ierr = MatAssemblyBegin( C, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
    ierr = MatAssemblyEnd( C, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
    ierr = MatDuplicate( C, MAT_COPY_VALUES, &newC );  CHKERRQ(ierr);
    if( scalar != 1.0 ) {
      ierr = MatScale( &scalar, newC );  CHKERRQ(ierr);
    }
  }
  else newC = NULL;

  data.setTypeName(matrix_name);
  PromCRMat *outmat = new PromCRMat(newA,newC);

  if( outmat == NULL ){ ERR_RETURN("setTypeName == NULL");}
  data.setDataPtr( (void*)outmat );

  return 0;
}

//PETSc_LinSysCore::sumInMatrix===============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::sumInMatrix"
int PETSc_LinSysCore::sumInMatrix( double scalar, const Data& data ) 
{
  PRINTD(scalar);
  int ierr;

  if( strcmp(matrix_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Matrix'.");
  }
  PromCRMat *AA = (PromCRMat*)data.getDataPtr(); 
  Mat A = K_->A_, C = K_->C_, inA = AA->A_, inC = AA->C_;
  if( (C==NULL) != (inC==NULL) ){ERR_RETURN("C==NULL != inC==NULL");}

  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
  ierr = MatAssemblyBegin( inA, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
  ierr = MatAssemblyEnd( inA, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
  //   Computes Y = a*X + Y. int MatAXPY(Scalar *a,Mat X,Mat Y)
  ierr = MatAXPY( &scalar, inA, A, SAME_NONZERO_PATTERN );  CHKERRQ(ierr);
  dirty_system_ = TRUE;

  if( C != NULL ){ 
    ierr = MatAXPY( &scalar, inC, C,  SAME_NONZERO_PATTERN );  CHKERRQ(ierr); 
    ierr = MatAssemblyBegin( C, MAT_FINAL_ASSEMBLY );  CHKERRQ(ierr);
    ierr = MatAssemblyEnd( C, MAT_FINAL_ASSEMBLY );   CHKERRQ(ierr);
  }

  return 0;
}

//PETSc_LinSysCore::getRHSVectorPtr==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getRHSVectorPtr"
int PETSc_LinSysCore::getRHSVectorPtr(Data& data) 
{
  PRINT;
  data.setTypeName( vector_name );
  data.setDataPtr( (void*)b_[currentRHS_] );

  return 0;
}

//PETSc_LinSysCore::copyInRHSVector===========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::copyInRHSVector"
int PETSc_LinSysCore::copyInRHSVector(double scalar, const Data& data) 
{
  PRINTD(scalar);
  int ierr;
  if( strcmp(vector_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Vector'.");
  }
  PromCRVec *yy = (PromCRVec*)data.getDataPtr(), *bb = b_[currentRHS_]; 

  ierr = VecCopy( yy->x_, bb->x_ );  CHKERRQ(ierr);
  if( bb->p_ != NULL ){
    ierr = VecCopy( yy->p_, bb->p_ );  CHKERRQ(ierr);
  }

  return 0;
}

//PETSc_LinSysCore::copyOutRHSVector==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::copyOutRHSVector"
int PETSc_LinSysCore::copyOutRHSVector( double scalar, Data& data ) 
{
  PRINTD(scalar);
  int ierr;
  PromCRVec *bb = b_[currentRHS_], *outvec; Vec x,p;
  
  ierr = VecDuplicate( bb->x_, &x );  CHKERRQ(ierr);
  ierr = VecCopy( bb->x_, x );  CHKERRQ(ierr);
  if( scalar != 1.0 ) {
    ierr = VecScale( &scalar, x );  CHKERRQ(ierr);
  }
  if( bb->p_ != NULL ) {
    ierr = VecDuplicate( bb->p_, &p );  CHKERRQ(ierr);
    ierr = VecCopy( bb->p_, p );  CHKERRQ(ierr);
    if( scalar != 1.0 ) {
      ierr = VecScale( &scalar, p );  CHKERRQ(ierr);
    }
  }
  else p = NULL;
  
  outvec = new PromCRVec(x,p);
  data.setTypeName( vector_name );
  data.setDataPtr( (void*)outvec );  

  return 0;
}

//PETSc_LinSysCore::sumInRHSVector============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::sumInRHSVector"
int PETSc_LinSysCore::sumInRHSVector( double scalar, const Data& data ) 
{
  PRINTD(scalar);
  int ierr;
  if( strcmp(vector_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Vector'.");
  }
  PromCRVec *xx = (PromCRVec*)data.getDataPtr(), *bb = b_[currentRHS_];
  
  // Computes y = alpha x + y. 
  ierr = VecAXPY( &scalar, xx->x_, bb->x_ );  CHKERRQ(ierr);
  if( xx->p_ != NULL ){ 
    ierr = VecAXPY( &scalar, xx->p_, bb->p_ );  CHKERRQ(ierr); 
  }

  return 0;
}

//PETSc_LinSysCore::destroyVectorData=========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::destroyVectorData"
int PETSc_LinSysCore::destroyVectorData(Data& data) 
{
  PRINT;
  if( strcmp(vector_name, data.getTypeName()) != 0 ){
    ERR_RETURN("data's type string not 'PETSc_Vector'.");
  }
  PromCRVec *xx = (PromCRVec*)data.getDataPtr();

  delete xx;
  data.setDataPtr( NULL );  

  return 0;
}

//PETSc_LinSysCore::setNumRHSVectors==========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setNumRHSVectors"
int PETSc_LinSysCore::setNumRHSVectors( int numRHSs, const int* rhsIDs ) 
{
  PRINT1(numRHSs); 
  int ierr, ii, oldnumrhs = numRHSs_;
  
  numRHSs_ = numRHSs;
  if ( numRHSs == oldnumrhs ) ; // noop
  else if( numRHSs == 0 ){ // deleting ???
    if( rhsIDs_ == NULL ){ ERR_RETURN("rhsIDs_ == NULL");} 
    PetscFree(rhsIDs_); rhsIDs_ = NULL; 
    for( ii = 0 ; ii < oldnumrhs ; ii++ ){
      PromCRVec *bb = b_[ii]; 
      if( bb == NULL ){ ERR_RETURN("bb == NULL");}
      delete bb; b_[ii] = NULL;
    }
    PetscFree( b_ ); b_ = NULL;
  }
  else if( oldnumrhs > numRHSs ) { // reducing 
    for( ii = numRHSs ; ii < oldnumrhs ; ii++ ){
      PromCRVec *bb = b_[ii]; b_[ii] = NULL; delete bb;
    }     
  }
  else { // increasing
    PetscFree( rhsIDs_ );
    ierr = PetscMalloc(sizeof(int)*numRHSs_,&rhsIDs_);CHKERRQ(ierr);
    PromCRVec** oldBarr = b_;
    ierr = PetscMalloc(numRHSs_*sizeof(PromCRVec*),&b_);CHKERRQ(ierr);
    for( ii = 0 ; ii < oldnumrhs ; ii++ ) b_[ii] = oldBarr[ii];
    PetscFree( oldBarr ); ii = oldnumrhs;
    do{ 
      Vec x, p = NULL;
      ierr = VecCreateMPI(FEI_PETSc_Comm_,numLocalRowsPrim_, 
			  numGlobalRowsPrim_,&x);
      CHKERRQ(ierr);  ierr = VecSet( &zero, x );  CHKERRQ(ierr);  
      if( numGlobalRowsLm_ ) { // LM Vecs
	ierr = VecCreateMPI(FEI_PETSc_Comm_,numLocalRowsLm_,numGlobalRowsLm_, 
			    &p);
	CHKERRQ(ierr);  ierr = VecSet( &zero, p );  CHKERRQ(ierr);
      }
      b_[ii] = new PromCRVec( x, p ); 
    }while( ++ii < numRHSs_ );
  }
  for( ii = 0; ii < numRHSs_; ii++ ) rhsIDs_[ii] = rhsIDs[ii];   

  return 0;
}

//PETSc_LinSysCore::setRHSID==================================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setRHSID"
int PETSc_LinSysCore::setRHSID( int rhsID ) 
{
  PRINT1(rhsID); 
  for(int i=0; i<numRHSs_; i++){
    if( rhsIDs_[i] == rhsID ){
      currentRHS_ = i;
      return 0;
    }
  }
  ERR_RETURN("rhsID not found.");
}

// PETSc_LinSysCore::formResidual ==========================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::formResidual"
int PETSc_LinSysCore::formResidual( double* values, int len )
{
  // explicit residual
  PRINT1(len);
  int ierr,ii,jj; PromCRVec *bb = b_[currentRHS_]; 
  double *pv;
  Vec xx, pp = NULL;
  if( bb == NULL ){ERR_RETURN("no RHS vector");}
  if( len != numLocalRows_ ){ ERR_RETURN("len != numLocalRows_"); }

  // rr =: Ax - b  
  ierr = VecDuplicate( x_->x_, &xx );  CHKERRQ(ierr);
  if( x_->p_  != NULL ) { ierr = VecDuplicate( x_->p_, &pp );  CHKERRQ(ierr); }
  PromCRVec *rr = new PromCRVec( xx, pp ); 
  ierr = K_->Mult( x_, rr ); CHKERRQ(ierr);
  ierr = rr->AXPY( &mone, bb ); CHKERRQ(ierr);
  // prim
  ierr = VecGetArray( rr->x_, &pv );    CHKERRQ(ierr);
  for( ii = 0 ; ii < numLocalRowsPrim_ ; ii++ ){
    values[ii] = -pv[ii]; // b - Ax = -(Ax - b)
  }
  ierr = VecRestoreArray( rr->x_, &pv ); CHKERRQ(ierr);
  // lm
  if( rr->p_ != NULL ) {
    ierr = VecGetArray( rr->p_, &pv );    CHKERRQ(ierr);
    for( jj = 0 ; jj < numLocalRowsLm_ ; jj++, ii++ ){
      values[ii] = -pv[jj]; // b - Ax = -(Ax - b)
    }
    ierr = VecRestoreArray( rr->p_, &pv ); CHKERRQ(ierr);
  }
  if( ii != numLocalRows_ ){ ERR_RETURN("ii != numLocalRows_"); }
  delete rr;

  return 0;
}

//PETSc_LinSysCore::getSolution===============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getSolution"
int PETSc_LinSysCore::getSolution( double* answers, int len)
{
  //
  //The caller must allocate the memory for 'answers' and eqnNumbers,
  //and len must be set to the right value -- i.e., len should equal
  //numLocalRows_.
  // 
  PRINT;
  if (len != numLocalRows_) {ERR_RETURN("len != numLocalRows_");}
  int ierr, jj, ii;
  double *vv; Vec xx = x_->x_, pp = x_->p_;
  if( len != numLocalRows_ ){ ERR_RETURN("len != numLocalRows_"); }

  ierr = VecGetArray( xx, &vv ); CHKERRQ(ierr);
  for( ii = 0 ; ii < numLocalRowsPrim_ ; ii++ ) {
    answers[ii] = vv[ii];
  }
  ierr = VecRestoreArray( xx, &vv );  CHKERRQ(ierr);
  // lm
  if( pp != NULL ) {
    ierr = VecGetArray( pp, &vv );    CHKERRQ(ierr);
    for( jj = 0 ; jj < numLocalRowsLm_ ; jj++, ii++ ){
      answers[ii] = vv[jj]; // b - Ax = -(Ax - b)
    }
    ierr = VecRestoreArray( pp, &vv ); CHKERRQ(ierr);
  }
  if( ii != numLocalRows_ ){ ERR_RETURN("ii != numLocalRows_"); }

  return 0;
}

//PETSc_LinSysCore::getSolnEntry==============================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getSolnEntry"
int PETSc_LinSysCore::getSolnEntry( int eqnNumber, double& answer ) 
{
  //
  //This function returns a single solution entry, the coefficient for
  //equation number eqnNumber.
  //
  //PRINT;
  if( localStartRow_ > eqnNumber || eqnNumber > localEndRow_ ){
    ERR_RETURN("eqnNumber out of range");
  }
  int ierr, myeq, lamb, proc;
  double *vv;
  ierr = eq_map_private( eqnNumber, myeq, lamb, proc ); CHKERRQ(ierr);
  if( proc!=thisProc_ ){ ERR_RETURN("pe != mype");}
  if( !lamb ) {
    ierr = VecGetArray( x_->x_, &vv ); CHKERRQ(ierr);
    answer = vv[myeq - localStartRowPrim_]; 
    ierr = VecRestoreArray( x_->x_, &vv ); CHKERRQ(ierr);
  }
  else {
    ierr = VecGetArray( x_->p_, &vv ); CHKERRQ(ierr);
    answer = vv[myeq - localStartRowLm_]; 
    ierr = VecRestoreArray( x_->p_, &vv ); CHKERRQ(ierr);
  }

  return 0;
}

#include <src/ksp/ksp/kspimpl.h>
static int setRegularization( Mat A, const Mat C, Vec xx, const double alpha,
			      Vec dd, const int startLmRow, const int numLmRow,
			      const int startPrimRow,const int numlocPrimRows);
static int AddRegularization( Mat A, const double alpha, const int ncols, 
			      int colids[], const double colv[] );

// PETSc_LinSysCore::launchSolver ==========================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::launchSolver"
int PETSc_LinSysCore::launchSolver( int& solveStatus, int& iterations ) 
{
  //
  //This function does any last-second setup required for the
  //linear solver, then goes ahead and launches the solver to get
  //the solution vector.
  //Also, if possible, the number of iterations that were performed
  //is stored in the iterations_ variable.
  //
  PRINT;
  int ierr, its, tt, outerits; 
  KSP sles;
  PromCRVec *bb = b_[currentRHS_];  assert(bb!=NULL);
  Mat A = K_->A_, C = K_->C_; 
  Vec u = x_->x_, f = bb->x_, gg = bb->p_, p = x_->p_;

  // finish BCs
  ierr = applyBCs_private(); CHKERRQ(ierr);

  // create
  ierr = KSPCreate( FEI_PETSc_Comm_, &sles ); CHKERRQ(ierr);

  // get solver options
  ierr = KSPSetFromOptions( sles ); CHKERRQ(ierr);

  // set matrix in KSP
  ierr = KSPSetOperators( sles, A, A, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);

  MPI_Allreduce( &dirty_guess_, &tt, 1, MPI_INT, MPI_MAX, FEI_PETSc_Comm_ );
  dirty_guess_ = tt;   // see if any have initial guess  
  if( C == NULL ) {
    // initial guess f := f - Ax -- p is done latter
    if( dirty_guess_ ) { 
      ierr = MatMult( A, init_guess_->x_, u );        CHKERRQ(ierr);
      ierr = VecAXPY( &mone, u, f );                  CHKERRQ(ierr);
    }
    // solve
    ierr = VecSet( &zero, u );  CHKERRQ(ierr);
    ierr = KSPSetRhs( sles, f );       CHKERRQ(ierr);
    ierr = KSPSetSolution( sles, u );         CHKERRQ(ierr);
    ierr = KSPSolve( sles ); solveStatus = ierr; 
    ierr = KSPGetIterationNumber( sles, &its ); CHKERRQ(ierr); 
    iterations = its;
    CHKERRQ(ierr); 
    // add in init. guess
    if( dirty_guess_ ) {
      ierr = VecAXPY( &one, init_guess_->x_, x_->x_ );  CHKERRQ(ierr);
    }
  }
  else { 
    // lagrange solve
    if( x_->p_ == NULL ){ ERR_RETURN("eqnNumber out of range"); }
    if( Pen_stiff_void_ == NULL ){ ERR_RETURN("Pen_stiff_ == NULL"); }
    Vec PenStiff = (Vec)Pen_stiff_void_;
    double       normC, temp, normR, normR_0;
    const double rtol_0 = sles->rtol;
    Vec rr, Cu_g, wrkp;
    ierr = VecDuplicate( u, &rr );  CHKERRQ(ierr);
    ierr = VecDuplicate( p, &wrkp );  CHKERRQ(ierr);
    ierr = VecDuplicate( p, &Cu_g );  CHKERRQ(ierr);
    // Get and set regularization
    if( dirty_system_ == TRUE ) {
      ierr = setRegularization( A, C, u, 5.0, PenStiff, localStartRowLm_,
				numLocalRowsLm_, localStartRowPrim_,
				numLocalRowsPrim_ );             CHKERRQ(ierr);
      dirty_system_ = FALSE;
      // assemble -- collective
      matrixLoadComplete(); 
    }
    // need to get current solution to start: lambda = Pen * C * x_real
    ierr = VecNorm( f, NORM_2, &normR_0 );               CHKERRQ(ierr);
    ierr = VecNorm( gg, NORM_2, &temp );  CHKERRQ(ierr);
    normR_0 = sqrt(normR_0*normR_0 + temp*temp);	// |b|_2    
    sles->rtol = 1.e-3; // varibale tol
    // r := f - (Au + C'p) -- vector used in primary Ax=f solves
    if( dirty_guess_ ){ 
      ierr = MatMult( A, init_guess_->x_, u );                CHKERRQ(ierr);
      ierr = MatMultTransposeAdd( C, init_guess_->p_, u, rr ); CHKERRQ(ierr);
      ierr = VecAYPX( &mone, f, rr );            CHKERRQ(ierr);
      ierr = VecCopy( init_guess_->p_, p );     CHKERRQ(ierr); 
      ierr = VecCopy( init_guess_->x_, u );     CHKERRQ(ierr);
    }
    else { 
      ierr = VecCopy( f, rr ); CHKERRQ(ierr); 
      ierr = VecSet( &zero, u );                           CHKERRQ(ierr);
      ierr = VecSet( &zero, p );                           CHKERRQ(ierr);
    }
    Vec uu = init_guess_->x_;
    for( outerits = its = 0; outerits < 10 ; outerits++ ) {
      // x <- A^-1 ( f - C' * lambda )
      ierr = VecSet( &zero, uu );                           CHKERRQ(ierr);
      ierr = KSPSetSolution( sles, uu ); CHKERRQ(ierr);
      ierr = KSPSetRhs( sles, rr );       CHKERRQ(ierr);
      ierr = KSPSolve( sles );         CHKERRQ(ierr);
      ierr = KSPGetIterationNumber( sles, &tt ); CHKERRQ(ierr); 
      ierr = VecAYPX( &one, uu, u );    CHKERRQ(ierr); // x_{i+1} := x_i + ...
      if( tt < 0 ) tt = -tt; its += tt;
      // increment lagrange mult. solution -- lam <- lam + Pen*(C*u - g)
      ierr = MatMult( C, u, Cu_g );       CHKERRQ(ierr);
      ierr = VecAXPY( &mone, gg, Cu_g );   CHKERRQ(ierr);
      ierr = VecPointwiseMult( Cu_g, PenStiff, wrkp );CHKERRQ(ierr);
      ierr = VecAYPX( &one, wrkp, p );    CHKERRQ(ierr);
      
      // residual
      ierr = MatMult( A, u, uu );   CHKERRQ(ierr);	
      ierr = MatMultTransposeAdd( C, p, uu, rr); CHKERRQ(ierr); // (Au + C'p)
      ierr = VecAYPX( &mone, f, rr );  CHKERRQ(ierr);  // r := f - (Au + C'p)

      // test convergence
      ierr = VecNorm( Cu_g, NORM_2, &normC );             CHKERRQ(ierr);      
      if( (temp = normC * 1.e0) < sles->rtol ) sles->rtol=temp; // new tolerence
      PetscPrintf(FEI_PETSc_Comm_,
		  "\t[%d]%s %d) %d iters,|Cu_g|=%9.2e rtol=%9.2e",
		  thisProc_, __FUNC__, outerits, tt, normC,  sles->rtol );
      if( normC < rtol_0 ||1){ // test #1 
	ierr = VecNorm( rr, NORM_2, &normR );  CHKERRQ(ierr);	// |f|_2
	temp = sqrt(normR*normR + normC*normC);
	PetscPrintf(FEI_PETSc_Comm_,"[%d] -- |b-K*x|=%9.2e(%9.2e <? %9.2e)\n",
		    thisProc_, normR, temp, rtol_0 * normR_0 );
	if( temp < rtol_0 * normR_0 ){ its++; break; } // test convergence #2
      }  
      else PetscPrintf( FEI_PETSc_Comm_, "\n" );
    }
    solveStatus = 0; // (temp >= rtol_0 * normR_0); // coverged?
    iterations = its;
    // clean up
    sles->rtol = rtol_0;
    ierr = VecDestroy( rr ); CHKERRQ( ierr );
    ierr = VecDestroy( wrkp ); CHKERRQ( ierr );
    ierr = VecDestroy( Cu_g ); CHKERRQ( ierr );
    
  }
  // clean up
  ierr = KSPDestroy(sles); CHKERRQ(ierr);
  ierr = VecSet( &zero, init_guess_->x_ );  CHKERRQ(ierr);      
  dirty_guess_ = FALSE;
    
  if( verbose_ > 2 ) {
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(f,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  return 0;
}


//PETSc_LinSysCore::eq_map_private=====================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::eq_map_private"
int PETSc_LinSysCore::eq_map_private( const int geq, // FEI-LSC eq number 
				      int &myeq,     // my eq. number (out)
				      int &lamb,     // 
				      int &proc ) const 
{
  if( verbose_ > 2 ) PRINT1(geq); // this is too noisy
  int ierr, nPrim,prim0,eq0,lm0,localEq; 

  ierr = getProc_private( geq, proc );  CHKERRQ(ierr);
  if( proc_lmEq_ == NULL ){ lamb = 0; myeq = geq; } // no LMs
  else {
    eq0 = proc_globalEq_[proc];
    lm0 = proc_lmEq_[proc];
    prim0 = proc_primEq_[proc];
    nPrim = proc_primEq_[proc+1] - prim0; 
    // map
    localEq = geq - eq0;
    if( localEq < nPrim ){ // Primary eq 
      lamb = 0; 
      myeq = prim0 + localEq; 
    }
    else{ // LM
      lamb = 1;
      myeq = lm0 + (localEq - nPrim); // stsrtLm + localLm
    }
  }
  
  return 0;
}

//PETSc_LinSysCore::getProc_private=====================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::getProc_private"
int PETSc_LinSysCore::getProc_private( const int geq, int &proc )const
{ 
  if( geq < 0 || geq >= numGlobalRows_ ){
    PRINT3(geq,numProcs_,numGlobalRows_);
    ERR_RETURN("gid >=0 && gid < proc_gnode[numprocs]"); 
  }
  if(proc_globalEq_[numProcs_] != numGlobalRows_){ERR_RETURN("ng!=pg[npp]");}
  float fproc = (float)geq * (float)numProcs_ / (float)numGlobalRows_ + 0.5;
  proc = (int)fproc; 
  if( proc > numProcs_ ) proc = numProcs_;
  while( geq < proc_globalEq_[proc] || geq >= proc_globalEq_[proc+1] ) {
    if(  geq < proc_globalEq_[proc] ) proc--;
    else proc++;
  }
  if( proc < 0 || proc >= numProcs_ ){
    ERR_RETURN("proc < 0 || proc >= numprocs"); 
  }

  return 0;
}

// PETSc_LinSysCore::setRegularization =============================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::setRegularization"
int setRegularization( Mat A, const Mat C, Vec xx, const double alpha, Vec dd,
		       const int startLmRow, const int numLmRow,
		       const int startPrimRow, const int numlocPrimRows )
{ 
  int ierr, ii, jj, primeq, loceq, lmeq, ncols, *colids, N; 
  double temp, *vv, *diagv = NULL, norm; 
  Vec diag;
  ierr = MatGetOwnershipRange(A, &ii, &jj); CHKERRQ(ierr); 
  ierr =  MatGetSize( A, &ii, &jj );  CHKERRQ(ierr);       N = jj;
  //PRINT2( (int)alpha, N );

  ierr = VecDuplicate( xx, &diag );  CHKERRQ(ierr);
  ierr = MatGetDiagonal( A, diag );      CHKERRQ(ierr);
  ierr = VecNorm( diag, NORM_1, &norm ); CHKERRQ(ierr);  norm = norm/(double)N;
  ierr = VecGetArray( diag, &diagv );    CHKERRQ(ierr);
  for( ii = 0, lmeq = startLmRow ; ii < numLmRow ; ii++, lmeq++ ){
    ierr = MatGetRow( C, lmeq, &ncols, &colids, &vv );   CHKERRQ(ierr);
    for( jj = 0, temp = 0.0 ; jj < ncols ; jj++ ) {
      primeq = colids[jj];
      loceq = primeq - startPrimRow;
      if( loceq >= 0 && loceq < numlocPrimRows ) temp += diagv[loceq];
      else { temp = -1.0; break; }
    }
    if( temp == 0.0 ) { std::cerr<<"temp == 0.0 ?????"<<std::endl;return(-1); }
    else if ( temp == -1.0 ) temp = norm;
    else temp = alpha*(temp/(double)ncols); // average of diagonal K
    // make panalty "matrix"
    ierr = VecSetValues( dd, 1, &lmeq, &temp, INSERT_VALUES );  CHKERRQ(ierr);
    // regularize A += temp * vv * vv'
    ierr = AddRegularization( A, temp, ncols, colids, vv );     CHKERRQ(ierr);
    ierr = MatRestoreRow( C, lmeq, &ncols, &colids, &vv ); CHKERRQ(ierr);
  } 
  ierr = VecRestoreArray( diag, &diagv );    CHKERRQ(ierr);
  ierr = VecDestroy( diag ); CHKERRQ( ierr );
  
  return 0;
}

// PETSc_LinSysCore::AddRegularization ==============================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::AddRegularization"
int AddRegularization( Mat A, const double alpha, const int ncols,int colids[],
		       const double colv[] )
{
  //PRINT;
  int ierr, ii, jj = ncols*ncols;  
  double *arr,*vp; 
  ierr = PetscMalloc(sizeof(double)*jj,&vp);CHKERRQ(ierr);
  arr = vp;

  for( ii = 0 ; ii < ncols ; ii++ ) {
    double ti = alpha * colv[ii];
    for( jj = 0 ; jj < ncols ; jj++ ) *vp++ = ti * colv[jj];
  }
  ierr = MatSetValues( A, ncols, colids, ncols, colids, arr, ADD_VALUES );
  CHKERRQ(ierr);
  PetscFree(arr);

  return 0;
}


/***********************************************************************
 *                              PromCRVec
 **********************************************************************/
// PromCRVec::Norm() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRVec::Norm"
int PromCRVec::Norm( NormType type, PetscReal *val )
{
  PetscReal temp; int ierr;

  ierr = VecNorm( x_, type, val );                           CHKERRQ(ierr);
  if( p_ != NULL ){ ierr = VecNorm( p_, type, &temp );      CHKERRQ(ierr);}
  else temp = 0.0;
  
  if( type == NORM_INFINITY ) *val = (*val>temp) ? *val : temp;
  else if( type == NORM_1 ) *val += temp;
  else{ // NORM_2 
    *val = sqrt( *val * *val + temp * temp );
  }
  return 0;
}

static char seterrq_str2[] = "!p_";

// PromCRVec::AYPX() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRVec::AYPX"
int PromCRVec::AYPX( const double *alpha, PromCRVec *x )
{
  int ierr;

  ierr = VecAYPX( alpha, x->x_, x_ );                      CHKERRQ(ierr);
  if( x->p_ != NULL ) {
    if( p_ == NULL ) SETERRQ(1,seterrq_str2); 
    ierr = VecAYPX( alpha, x->p_, p_ );                      CHKERRQ(ierr);
  }
  return 0;
}

// PromCRVec::AXPY() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRVec::AXPY"
int PromCRVec::AXPY( const double *alpha, PromCRVec *x )
{
  int ierr;

  ierr = VecAXPY( alpha, x->x_, x_ );                      CHKERRQ(ierr);
  if( x->p_ != NULL ) {
    if( p_ == NULL ) SETERRQ(1,seterrq_str2); 
    ierr = VecAXPY( alpha, x->p_, p_ );                      CHKERRQ(ierr);
  }
  return 0;
} 

static char seterrq_str3[] = "x_ == NULL";

// PromCRVec::Set() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRVec::Set"
int PromCRVec::Set( const double alpha )
{
  int ierr;

  if( x_ == NULL ) SETERRQ(1,seterrq_str3);
  ierr = VecSet( &alpha, x_ ); CHKERRQ(ierr);
  if( p_ != NULL ) { ierr = VecSet( &alpha, p_ ); CHKERRQ(ierr); }
  
  return 0;
}

static char seterrq_str4[] = "!x->p_";
static char seterrq_str5[] = "!y->p_";

/***********************************************************************
 *                           PromCRMat
 **********************************************************************/
// PromCRMat::Mult() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRMat::Mult"
int PromCRMat::Mult( PromCRVec *x, PromCRVec *y )
{
  int ierr;

  if( C_ != NULL ) {
    if( x->p_ == NULL ) SETERRQ(1,seterrq_str4); 
    if( y->p_ == NULL ) SETERRQ(1,seterrq_str5); 

    Vec workx; 
    ierr = VecDuplicate( x->x_, &workx );  CHKERRQ(ierr);
    ierr = MatMult( A_, x->x_, workx );                 CHKERRQ(ierr); 
    // v3 = v2 + A' * v1.    
    ierr = MatMultTransposeAdd( C_, x->p_, workx, y->x_ ); CHKERRQ(ierr);
    ierr = MatMult( C_, x->x_, y->p_ );                     CHKERRQ(ierr);
    ierr = VecDestroy( workx ); CHKERRQ( ierr );
  }
  else {
    ierr = MatMult( A_, x->x_, y->x_ );                   CHKERRQ(ierr);
  }
  return 0;
}

// PromCRMat::ZeroEntries() ======================================
#undef __FUNC__    
#define __FUNC__ "PromCRMat::ZeroEntries"
int PromCRMat::ZeroEntries()
{
  int ierr;

  ierr = MatZeroEntries( A_ ); CHKERRQ(ierr); 
  if( C_ != NULL ) { ierr = MatZeroEntries( C_ ); CHKERRQ(ierr); }
  return 0;
}

/*** NOT implemetnted ********************************************************/

//PETSc_LinSysCore::setLoadVectors===========================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setLoadVectors"
int PETSc_LinSysCore::setLoadVectors(GlobalID elemBlock,
				     int numElems,
				     const GlobalID* elemIDs,
				     const double *const *load,
				     int numEqnsPerElem,
				     const int *const * eqnIndices)
{
  return 0;
}

//PETSc_LinSysCore::setStiffnessMatrices =====================================
#undef __FUNC__    
#define __FUNC__ "PETSc_LinSysCore::setStiffnessMatrices"
int PETSc_LinSysCore::setStiffnessMatrices(GlobalID elemBlock,
					   int numElems,
					   const GlobalID* elemIDs,
					   const double *const *const *stiff,
					   int numEqnsPerElem,
					   const int *const * eqnIndices)
{
  return 0;
}


// PETSc_LinSysCore::writeSystem ==========================================
#undef __FUNC__
#define __FUNC__ "PETSc_LinSysCore::writeSystem"
int PETSc_LinSysCore::writeSystem(const char* name)
{
  // easy to do with PETSc
  return 0;
}
