/*
  This is the base class for all load-balancer

  jaideep ray, 08/23/02
*/

#ifndef BaseLBHSeen
#define BaseLBHSeen

// CCA specific
#include "cca.h"
#include "mpi.h"

// abstract class for all load-balancers
#include "LoadBalancer.h"

// things i use
#include "EntityList.h"
#include "EntityListImpl.h"
#include "PartitionerFactory.h"

// Zoltan specific
#include "include/zoltan.h"

//C++ io
#include <iostream>

// C++ stl
#include <map>
#include <string>
using namespace std;

namespace ZoltanSpace
{
  class BaseLB : public virtual ::LoadPartitionerSpace::LoadBalancer
  {
    public :

    BaseLB(PartitionerFactory_JR *q, int index, MPI_Comm *A)
    {      
      pFac = q ; my_index = index ; ref = 1 ; incoming = outgoing = 0 ;
      my_zz = Zoltan_Create( *A );

      MPI_Comm_rank( *A, &rank) ;
      MPI_Allreduce( &rank, &master, 1, MPI_INT, MPI_MIN, *A) ;

      int nprocs ;
      MPI_Comm_size( *A, &nprocs) ;
      my_proc_list = new int [nprocs] ;

      Zoltan_Set_Param(my_zz, "RETURN_LISTS", "ALL");
      Zoltan_Set_Param(my_zz, "IMBALANCE_TOL", "1.1");
      Zoltan_Set_Param(my_zz, "AUTO_MIGRATE", "FALSE");

      is_init = false ;
    }

    virtual ~BaseLB()
    { 
      pFac->removeLB( my_index ) ;
      if (incoming != 0) delete incoming ;
      if (outgoing != 0) delete outgoing ;
      Zoltan_Destroy( &my_zz ) ;
      if ( my_proc_list != 0 ) delete [] my_proc_list ;
    }
    
    virtual int CreateBalancedPartition( bool *did_it_change, 
		        ::LoadPartitionerSpace::EntityList **ElementsToImport,
			::LoadPartitionerSpace::EntityList **ElementsToExport);
    
    virtual int EvaluateDecomposition(bool print_stats, int *nobjs_on_proc, 
				      float *obj_wt_on_proc, int *ncuts_on_proc,
				      float *cut_wgts_on_proc, int *nbndry, 
				      int *n_my_adjacent_procs);
    
    virtual int cleanup();
    
    Zoltan_Struct * get_zz() { return ( my_zz) ; }

    virtual void addRef() { ref++ ; } 

    virtual void deleteRef() { ref-- ; if (ref == 0) delete (this) ; }

    virtual int GetIndex() { return (my_index) ; }

    // These methods are to be implemented by specific load-balancers.

    virtual int IncrementalAssignPoint(double *coords, int ndims, int *proc) = 0;
    
    virtual int IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, int *nprocs,
				     int **proc_list) = 0 ;
    
    virtual int GetLBInfo(char **name, int *index) = 0;
    
    virtual int SetParameter(char *key, char *val) = 0 ;
    virtual int GetParameter(char *key, char **val) = 0 ;

    virtual int SetParameter(char *key, double d) = 0 ;
    virtual int GetParameter(char *key, double *d) = 0 ;
    
    virtual int SetParameter(char *key, bool b) = 0 ;
    virtual int GetParameter(char *key, bool *b) = 0 ;
    
    virtual int SetParameter(char *key, int i) = 0 ;
    virtual int GetParameter(char *key, int *i) = 0;
    
    virtual int PrintKeywordsToScreen() = 0 ;

    protected :
      
    // Timer : Wall/CPU ; default  wall
    virtual int SetCommonParameter(char *key, char *val) ;
    virtual int GetCommonParameter(char *key, char **val) ;
    
    /* ImbalanceTolerance : a number; 1.2 is good (also default) */
    virtual int SetCommonParameter(char *key, double d) ;
    virtual int GetCommonParameter(char *key, double *d) ;

    /* Automigrate : hard-coded to false.
       Deterministic : false (default true)
       UseMachFile : hard-coded to false now
    */
    virtual int SetCommonParameter(char *key, bool b) ;
    virtual int GetCommonParameter(char *key, bool *b) ;

    /* NumGidEntries :  a number, default 1
       NumLidEntries :  a number, default 1
       ObjWtDim      :  a number {0 / 1}, default 0
       EdgeWtDim     :  a number {0 / 1}, default 0
       DebugLevel    :  a number > 0 ; default 1.
       DebugProc     :  a number > 0 ; default 0.
       CommWtDim     :  a number > 0; default 1
    */
    virtual int SetCommonParameter(char *key, int i) ;
    virtual int GetCommonParameter(char *key, int *i);

    // Print all keywords to screen
    virtual int PrintCommonKeywordsToScreen() ;

    int my_index, rank, master ;
    int *my_proc_list ;
    struct Zoltan_Struct *my_zz ;
    EntityListImpl *incoming, *outgoing ;
    PartitionerFactory_JR *pFac ;

    private :

    bool is_init ;
    void init() ;
    int ref ;

    map<string, string> props ;
    map<string, string> translation ;
  };
};
#endif

