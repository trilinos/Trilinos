#include "OctTree.h"
#include "PartitionerFactory.h"

void ZoltanSpace::OctTree_LB::init()
{
  prop["OCT_DIM"] = "3" ;
  prop["OCT_METHOD"] = "2";
  prop["OCT_MINOBJECTS"] = "10" ;
  prop["OCT_MAXOBJECTS"] = "40" ;
  prop["OCT_OUTPUT_LEVEL"] = "0" ;

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

  is_init = true ;
}

int ZoltanSpace::OctTree_LB::IncrementalAssignPoint(double *coords, int ndims, int *proc)
{
  cerr << "  ZoltanSpace::OctTree_LB::IncrementalAssignPoint() "
       << " Not implemented for " << myname << " partitioners. Return with -1"
       << endl ;
  return(-1) ;
}
    
int ZoltanSpace::OctTree_LB::IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, 
					      int *nprocs, int **proc_list)
{
  cerr << "  ZoltanSpace::OctTree_LB::IncrementalAssignBox() "
       << " Not implemented for " << myname << " partitioners. Return with -1"
       << endl ;
  return(-1) ;
}
