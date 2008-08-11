#include "ParMetis.h"
#include "PartitionerFactory.h"

void ZoltanSpace::ParMetis_LB::init()
{
  prop["PARMETIS_METHOD"] = "RepartGDiffusion" ;
  prop["PARMETIS_OUTPUT_LEVEL"] = "0";
  prop["PARMETIS_COARSE_ALG"] = "2" ;
  prop["PARMETIS_SEED"] = "15" ;
  prop["PARMETIS_ITR"] = "100" ;
  prop["PARMETIS_USE_OBJ_SIZE" ] = "1";
  prop["CHECK_GRAPH" ] = "1";
  prop["SCATTER_GRAPH" ] = "1";
  
  ::map<string,string>::iterator p;
  for(p = prop.begin() ; p != prop.end(); p++)
  {
    char *c1 = const_cast<char *>( ((*p).first).c_str() );
    char *c2 = const_cast<char *>( ((*p).second).c_str() );
    int res = Zoltan_Set_Param(BaseLB::my_zz, c1, c2) ;
  }

  // Set up all the query functions that this load balancer needs. Invoke
  // the right methods on the PartitionerFactory - that has a list of 
  // Ports etc.

  BaseLB::pFac->SetZoltanQueriesEntityInfoPort( BaseLB::my_zz ) ;
  BaseLB::pFac->SetZoltanQueriesGeomInfoPort( BaseLB::my_zz ) ;
  BaseLB::pFac->SetZoltanQueriesEdgeInfoPort( BaseLB::my_zz ) ;

  is_init = true ;
}

int ZoltanSpace::ParMetis_LB::IncrementalAssignPoint(double *coords, int ndims, int *proc)
{
  cerr << "  ZoltanSpace::ParMetis_LB::IncrementalAssignPoint() "
       << " Not implemented for " << myname << " partitioners. Return with -1"
       << endl ;
  return(-1) ;
}
    
int ZoltanSpace::ParMetis_LB::IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, 
					      int *nprocs, int **proc_list)
{
  cerr << "  ZoltanSpace::ParMetis_LB::IncrementalAssignBox() "
       << " Not implemented for " << myname << " partitioners. Return with -1"
       << endl ;
  return(-1) ;
}
