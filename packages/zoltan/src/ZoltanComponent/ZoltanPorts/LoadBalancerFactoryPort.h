/*
  A Port that fronts for a component that makes LoadBalancers
  Hasn't been debated on.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef LoadBalancerFactoryPortSeen
#define  LoadBalancerFactoryPortSeen

#include "cca.h"
#include "LoadBalancer.h"
#include "mpi.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace
{
  class LoadBalancerFactory : public virtual CONV_NS(Port)
  {

    public :

    LoadBalancerFactory()  : CONV_NS(Port)() {} 

    virtual ~LoadBalancerFactory() {} 

    /// Create a partitioner-virtual parallel machine pair
    /** create a partitioner and associate it with a parallel communicator
     *  @param A : pointer to a MPI communicator
     *  @param name : what kind of a partitioner do you want
     *  @return -1 : don't know about the partitioner you want
     *          -2 : NULL communicator passed.
     */
    virtual LoadBalancer *CreateLB_VM_pair(MPI_Comm *A, char *name) = 0 ;

    /// Destroy the load-balancer; don't need it any more.
    virtual void Destroy_LB_VM_pair( LoadBalancer *p)  = 0 ;
  } ;
};
#endif
