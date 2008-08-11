/*
  This is the implementation of the RCB wrapper.

  jaideep ray, 08/26/02
*/

#ifndef RCB_LBHSeen
#define RCB_LBHSeen

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
  class RCB_LB : public virtual BaseLB
  {
    public :

    RCB_LB(PartitionerFactory_JR *q, int index, MPI_Comm *A) : BaseLB(q, index, A)
    {      
      myname = "RCB" ; is_init = false ;
      Zoltan_Set_Param(BaseLB::my_zz, "LB_METHOD", "RCB");

      init() ;
    }

    virtual ~RCB_LB() { is_init = false ; }
    
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
	      2.   RCB_OVERALLOC : Extra temporary memory allocated by RCB to hold
	           objects when object-to-processor assignments change. 1.0 = 
		   no extra memory, 1.5 - 50 % more.  By default, no extra memory.
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
	      8. RCB_REUSE : flag to indicate whether to use the previous invocation's
	         cuts to get an initial guess for this invocation. 0 = don't use previous
		 cuts, 1 = use them. Default : 0
	      9. RCB_OUTPUT_LEVEL : flag controlling the amount of timing and diagnostic
	         info required. 0 = no output, 1 = print a summary, 2 = print report 
		 for each proc. Default : 0
	     10. CHECK_GEOM : flag controlling error-checking for inputs and outputs.
	         0 = no checking, 1 = check it. Default : 0
	     11. KEEP_CUTS : flag controlling whether to preserve the current 
	         decomposition. Helps if you need to add in an object later in a balanced
		 fashion. 0 = do not preserve the decomposition, 1 = do so. Default : 0.
	     12. RCB_LOCK_DIRECTIONS : flag to determine if the direction of cuts are to
	         preserved after 1st invocation, so that they can be repeated over-and-over
		 again in subsequent invocations.
	     13. RCB_SET_DIRECTIONS : if this flag is set, all the cuts in one direction
	         are done together. The cuts are determined first and then the value to
		 this parameter is checked regarding how to proceed. 0 = don't order
		 cuts, 1 = xyz, 2=xzy, 3=yzx, 4=yxz, 5=zxy, 6=zyx. Default : 0
             14. RCB_RECTILINEAR_BLOCKS :  flag controlling the shape of the resulting
	         region. If specified, when a cut is made, dots located on the cut are
		 moved to the same side of the cut. the region is then rectilinear,
		 but the load-balance is less uniform. Alternatively, the points can
		 be moved individually. 0 = move dots individually 1 = move as a goup.
		 Default : 0
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

