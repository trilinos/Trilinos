/*
  An abstract interface to allow the configuring of a partitioner.
  Hasn't been debated on.
  I am deriving this from a Port just in case someone wants to have a 
  1-load-balancer component i.e a RSB component.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef LoadBalancerInterfaceHeaderSeen
#define  LoadBalancerInterfaceHeaderSeen

#include "cca.h"
#include "EntityList.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace 
{
  class LoadBalancer : public virtual CONV_NS(Port)
   {

     public :

     LoadBalancer() : CONV_NS(Port)() {}

     virtual ~LoadBalancer() {}

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
     
     virtual int SetParameter(char *key, char *val) = 0 ;
     virtual int GetParameter(char *key, char **val) = 0 ;

     /** To set/get what the imbalance tolerance should be 
	 In case of GetParameter(), you allocate a double and I'll
	 fill it in.
	 @param key = ImbalanceTolerance (case-sensitive!)
	 @param val = a double precision number specifying the imbalance.
	 Each load-balancer defines its own imbalance metric but roughly
	 should be max_load / av_load, where these numbers are
	 calculated for a given mesh distribution across all processors.
	 E.g. ImbalanceTolerance = 1.2 is nice.
	 @return - 1 if the keyword is not found - perhaps the load-balancer
	 doesn't allow the setting of the imbalance tolerance
     */
     
     virtual int SetParameter(char *key, double d) = 0 ;
     virtual int GetParameter(char *key, double *d) = 0 ;

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
     
     virtual int SetParameter(char *key, bool b) = 0 ;
     virtual int GetParameter(char *key, bool *b) = 0 ;

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
	 @param value = an integer number
	 @return -1 if the key is not offered by the load-balancer for modification.
     */

     virtual int SetParameter(char *key, int i) = 0 ;
     virtual int GetParameter(char *key, int *i) = 0;

     /** Print to screen the list of keywords.
	 @return 0 if no error.
     */
     virtual int PrintKeywordsToScreen() = 0 ;
     //@}

     /// Methods that suggest a better partition of the current mesh
     /** These methods create a new partition and return it to you. You
	 may ignore this suggestion without any problems at all.
     */
     //@{
     /// Partition a mesh
     /** Pass in the current mesh and get a better partitioned mesh out. You are
	 welcome to ignore this new partition. The load-balancer will not actually
	 migrate the data for you.
	 @param did_it_change : Is the new partition different from the previous one ?
	        you allocate a single bool and pass its pointer in; I fill it.
	 @param ListToImport : list of elements to import from other procs.
	 @param ListToExport : list of my elements to export to other procs
	 @return 0 if successful, negative if error.
     */
     virtual int CreateBalancedPartition( bool *did_it_change, 
					  EntityList **ElementsToImport,
					  EntityList **ElementsToExport) = 0 ;

     /// Evaluate a partition.
     /** One you have called CreateBalancedPartition, you can immediately evaluate
	 it. Don't wait since someone might cleanup and remove all memory 
	 @param print_stats : Print statiistics to stdout
	 @param nobjs_on_proc : You allocate an int; I fill it with the number of
	        entities (mesh points/cells etc) on this processor
	 @param obj_wt_on_proc :  Total computational load for this proc
	 @param ncuts_on_proc : You allocate an int, I fill it with the number of
	        entity-to-entity connectivities I have had to cut to create this 
		partition.
	 @param cut_wgths_on_proc : A measure of the communications cost during
	        ghost-cell updates.
	 @param nbndy : You alloctate an int, I fill it with the number of 
	                domain boundary objects on this proc.
	 @param n_my_adjacent_procs : You allocate an int, I fill with the 
	        number of adjacent procs this processor need to communicate with.
        @return -1 if something goes wrong.
     */
     virtual int EvaluateDecomposition(bool print_stats, int *nobjs_on_proc, 
				       float *obj_wt_on_proc, int *ncuts_on_proc,
				       float *cut_wgts_on_proc, int *nbndry, 
				       int *n_my_adjacent_procs) = 0 ;
     
     /// Done with this partition, remove all memory.
     /** Done with this partition, remove all memeory and get ready for the
	 next partitioning task.
     */
     virtual int cleanup() = 0 ;
     //@}

     /// Incremental assignment methods
     /** If an extra point or area had to be added to the partitioned domain,
	 where would it go ?
     */
     //@{
     /// If an additional point is to be added
     /** If an additional point needs to be added to the decomposed mesh,
	 which processor would it go to ?
	 @param coords : pointer to an array you allocated, containing the
	 coordinates of the new point
	 @param ndims : length of the coords array
	 @param proc_no : pointer to the int you allocated which I will
	 fill with the processor this point should be on.
	 @return -1 : something went wrong.
     */
     virtual int IncrementalAssignPoint(double *coords, int ndims, int *proc) = 0 ;

     /// If a rectangular region needs to be added ....
     /** If a rectangular region needs to be added, which processors should it
	 go onto ?
	 @param lbbc : pointer to an array containing the lower bounding box corner
	               of the box.
	 @param ubbc : pointer to an array containing the upper bounding box corner
	               of the box
	 @param ndims : length of the lbbc & ubbc arrays
	 @param nprocs :  you allocate an int, I fill it with the number of
	                  processors the box will get distributed onto
	 @param proc_list : a pointer to my array containing which processors
	                    the box got distributed onto.
	 @return -1 : Error
     */
     virtual int IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, 
				      int *nprocs, int **proc_list) = 0 ;
     //@}


     /// Get load balancer info and index
     /** Gets the name and index of load-balancer
      *  @param  name : a pointer to my char array containing my name
      *  @param  index : a number that I carry to denote I am the index-th 
      *  load-balancer, in case I am in an environment with multiple load-balancers.
      *  @return -1 in case of error.
      */
     virtual int GetLBInfo(char **name, int *index) = 0 ;

     virtual void addRef() = 0 ;

     virtual void deleteRef() = 0 ;
   } ;
} ;
#endif

     



