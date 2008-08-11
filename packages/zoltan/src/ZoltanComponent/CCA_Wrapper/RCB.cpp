#include "RCB.h"
#include "PartitionerFactory.h"

void ZoltanSpace::RCB_LB::init()
{
  prop["RCB_OVERALLOC"] = "1.0" ;
  prop["RCB_REUSE"] = "0" ;
  prop["RCB_OUTPUT_LEVEL"] = "0";
  prop["CHECK_GEOM"] = "1" ;
  prop["KEEP_CUTS"] = "0" ;
  prop["RCB_LOCK_DIRECTIONS"] = "0";
  prop["RCB_SET_DIRECTIONS"] = "0" ;
  prop["RCB_RECTILINEAR_BLOCKS"] = "0" ;

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

int ZoltanSpace::RCB_LB::IncrementalAssignPoint(double *coords, int ndims, int *proc)
{
  int res = Zoltan_LB_Point_Assign(BaseLB::my_zz, coords, proc) ;
  
  if ( res == ZOLTAN_FATAL || res == ZOLTAN_MEMERR ) return(-1) ;

  return(0) ;
}
    
int ZoltanSpace::RCB_LB::IncrementalAssignBox(double *lbbc, double *ubbc, int ndim, 
					      int *nprocs, int **proc_list)
{
  double xmin, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0 ;
  xmin = lbbc[0] ; if (ndim > 1) ymin = lbbc[1] ; if (ndim > 2) zmin = lbbc[2] ;
  xmax = ubbc[0] ; if (ndim > 1) ymax = ubbc[1] ; if (ndim > 2) zmax = ubbc[2] ;

  int res = Zoltan_LB_Box_Assign(BaseLB::my_zz, xmin, ymin, zmin, xmax, ymax,
				 zmax, BaseLB::my_proc_list, nprocs) ;

  if ( res == ZOLTAN_FATAL || res == ZOLTAN_MEMERR ) return (-1) ;

  *proc_list = my_proc_list ;

  return(0) ;
}
