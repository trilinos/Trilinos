/*
  This is the implementation of the ParMetis wrapper.

  jaideep ray, 08/26/02
*/

#ifndef ParMetis_LBHSeen
#define ParMetis_LBHSeen

// CCA specific
#include "cca.h"
#include "mpi.h"

// abstract class for all load-balancers
#include "LoadBalancer.h"

// things i use
#include "PartitionerFactory.h"

// The base class I extend
#include "BaseLB.h"

// Zoltan specific
#include "include/zoltan.h"

// C++
#include <iostream>
#include <string>
#include <map>
using namespace std;

namespace ZoltanSpace
{
  class ParMetis_LB : public virtual BaseLB
  {
    public :

    ParMetis_LB(PartitionerFactory_JR *q, int index, MPI_Comm *A) : BaseLB(q, index, A)
    {      
      myname = "ParMetis" ; is_init = false ;
      Zoltan_Set_Param(BaseLB::my_zz, "LB_METHOD", "PARMETIS");

      init() ;
    }

    virtual ~ParMetis_LB() { is_init = false ; }
    
     /// Set/Get methods to configure a load balancer
     /** @name Load-balancer configuring methods 
      * These methods configure a load-balancer based on key-value
      * pairs. The actual key-value pairs are below. Note : For the
      * Get
      */

     //@{
     /**
	To set/get whether  to use CPU or Wall clock as timer.
	@param key = Timer
	@param value = Wall / CPU (case-sensitive !)
	@return -1 if the key was not found (and hence this load-balancer 
	does not allow the setting of a timer (to measure its own performance 
	time).
     */
    virtual int SetParameter(char *key, char *val) ;
    virtual int GetParameter(char *key, char **val) ;

     /** To set/get what the imbalance tolerance should be 
	 In case of GetParameter(), you allocate a double and I'll
	 fill it in.
	 @param key :
	      1.  ImbalanceTolerance (case-sensitive!) Each load-balancer 
	           defines its own imbalance metric but roughly should be 
		   max_load / av_load, where these numbers are calculated for 
		   a given mesh distribution across all processors.
		   E.g. ImbalanceTolerance = 1.2 is nice.
	 @param val = a double precision number. specifying the imbalance.
	 @return - 1 if the keyword is not found - perhaps the load-balancer
	 doesn't allow the setting of the imbalance tolerance
     */
    virtual int SetParameter(char *key, double d) ;
    virtual int GetParameter(char *key, double *d) ;

    /** To set/get whether I should (a) auto-migrate i.e move "elements"
	 for you to accomplish a load-balance (as if the load-balancer knows
	 anything about the mesh elements) (b) deterministic (i.e. gives the
	 same result everytime it runs - bit of a bother since asynchronous
	 communications cannot be done (c) if one should use a machine
	 description file for heterogeneous clusters. In case of GetParameter(),
	 you allocate a bool and I'll fill it in.
	 @param key = "AutoMigrate", "Deterministic, "UseMachFile" (case-sensitive)
	 @param val = true or false
	 @return -1 if this particular key is not offered for chnaging by user
     */
    virtual int SetParameter(char *key, bool b) ;
    virtual int GetParameter(char *key, bool *b) ;

     /** To set/get some very important parameter describing the mesh etc.
	 In case of GetParameter() method, you allocate an integer and I will
	 fill it in.
	 @param key : 
	      1. NumGidEntries : Entities in the mesh to be partitioned are
		 identified by unsigned int IDs, unique across all processors. 
		 These need not be single numbers; each "entity" (a mesh 
		 points/element) can be identified by a 1D array too, Q-elements 
		 long. NumGidEntried sets this integer value  Q.
	      2. NumLidEntries : Just as "entities" are identified by an ID globally,
	         they might be identified locally on a processors by another set of
		 unsigned ints (perhaps a direct mapping into a local-to-a-processor
		 array ?). This local ID can be a 1D array, R-elements long. Set "R"
		 using this keyword.
	      3. ObjWtDim : Each "entity" has a weight proportional to its CPU
	         requirements. Again, this can be a 1D array, S-elements long.
		 Set "S" using this keyword.
	      4. EdgeWtDim : Same as above, but Edge weights are proportional to
	         the communication requirements for a mesh point.
	      5. DebugLevel : Starting from 0, increasingly more debug info.
              6. DebugProc : Starting from 0, which proc dumps debug info. It will
	         dump info only for itself.
	      7. CommWtDim : length of an array which contains the communication costs
                 between procs.
	      8. PARMETIS_METHOD : The ParMetis method to be used - 9 available now
	         + PartKway - multilevel Kernighan-Lin partitioning
		 + PartGeom - space-filling curve (coordinates based)
		 + PartGeomKway - hybrid method, based on both PartKway and PartGeom
		   and needs graph and coordinates data.
		 + AdaptiveRepart - adaptive repartitioning, but only if ParMetis 3.0 
		   (or above) exists
		 + RepartLDiffusion - a local diffusion algorithm
		 + RepartGDiffusion - a global diffusion algorithm
		 + RepartRemap - multilevel partitioning with remap to minimize 
		   migration costs.
		 + RepartMLRemap - similar to above, but with additional multilevel
		   refinement.
		 + RefineKway - refine the current partitioning (balance)
		 Default : RepartGDiffusion
              9. PARMETIS_OUTPUT_LEVEL : amount of reporting to be done
	         0 = no report, 1 = timing info.  Default : 0
	     10. PARMETIS_COARSE_ALG : Coarse algorithm for PartKway. 1 = serial,
	         2 = parallel. Default : 2
	     11. PARMETIS_SEED : a random seed for ParMETIS. Default : 15
	     12. PARMETIS_ITR : Ratio of interprocessor communication time to 
	         redistribution time. A high value will emphasize reducing 
		 edge cuts, a low one will try to minimize changes between 2 
		 consecutive partitionings. A value of 100 -- 1000 is good. 
		 Used only for AdaptiveRepart. Default: 100.
	     13. PARMETIS_USE_OBJ_SIZE : Use the info re object sizes to estimate
	         migration costs. Only relevant for AdaptiveRepart. Default. 1
	     14. CHECK_GRAPH : Level or error checking for graph input. 0 = 
	         no checking, 1 = one-proc checking, 2 = global checking. 
		 The last is very slow; use only for debugging.
	     15. SCATTER_GRAPH : Scatter the graph data by distributing roughly
	         equal chunks of entities (graph vertices) to processors before
		 calling the partitioner. 0 = don't scatter, 1 = scatter only if
		 all entities are on a single-proc, 2 = scatter if there is at least
		 one proc with on entities, 3 = always scatter. Default = 1
	 @param value = an integer number
	 @return -1 if the key is not offered by the load-balancer for modification.
     */
    virtual int SetParameter(char *, int i) ;
    virtual int GetParameter(char *, int *i) ;

    /** Print to screen the list of keywords.
	@return 0 if no error.
     */
    virtual int PrintKeywordsToScreen() ;
    
    //@}

    virtual int IncrementalAssignPoint(double *coords, int ndims, int *proc) ;
    
    virtual int IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, 
				     int *nprocs, int **proc_list) ;
    
    virtual int GetLBInfo( char **name, int *index)
    {
      *name = const_cast< char * > (myname.c_str()) ;   *index = BaseLB::my_index;
    }

    private :

    void init() ;
    bool is_init ;

    string myname ;
    map<string, string> prop ;
  };
};
#endif

