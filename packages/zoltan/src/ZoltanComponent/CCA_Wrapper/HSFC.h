/*
  This is the implementation of the RCB wrapper.

  jaideep ray, 08/26/02
*/

#ifndef HSFC_LBHSeen
#define HSFC_LBHSeen

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
  class HSFC_LB : public virtual BaseLB
  {
    public :

    HSFC_LB(PartitionerFactory_JR *q, int index, MPI_Comm *A) : BaseLB(q, index, A)
    {      
      myname = "HSFC" ; is_init = false ;
      Zoltan_Set_Param(BaseLB::my_zz, "LB_METHOD", "HSFC");

      init() ;
    }

    virtual ~HSFC_LB() { is_init = false ; }
    
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
	     8. KEEP_CUTS : flag controlling whether to preserve the current 
	        decomposition. Helps if you need to add in an object later in a balanced
		fashion. 0 = do not preserve the decomposition, 1 = do so. Default : 0.
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

