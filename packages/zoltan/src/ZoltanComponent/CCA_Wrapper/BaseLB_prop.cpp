/*
  Implementation of setting of common properties for all load balancers.
  Jaideep Ray, SNL, 08/27/02
*/

#include "BaseLB.h"
#include <sstream>

// Timer : Wall/CPU ; default  wall
int ZoltanSpace::BaseLB::SetCommonParameter(char *key, char *val)
{

  if (is_init == false) init() ;

  string name(key), value(val) ;

  if ( translation.find(name) == translation.end() ) return(-1) ;

  string dummy = translation[ name ] ;
  props[ dummy ] = val ;

  return(0) ;
}

int ZoltanSpace::BaseLB::GetCommonParameter(char *key, char **val)
{

  if (is_init == false) init() ;

  string name(key);

  if ( translation.find(name) == translation.end() ) return(-1) ;

  string dummy = translation[ name ] ;
  *val = const_cast<char *> ((props[ dummy ]).c_str()) ;

  return(0) ;
}  
    
/* ImbalanceTolerance : a number; 1.2 is good (also default) */
int ZoltanSpace::BaseLB::SetCommonParameter(char *key, double d) 
{

  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;

  stringstream dummy ; dummy << d << ends ;
  string charform; dummy >> charform ;

  props[ ZoltanName ] = charform ;
  Zoltan_Set_Param( my_zz, const_cast<char *> (ZoltanName.c_str()), 
		    const_cast<char *>((props[ ZoltanName ]).c_str()) );

  return(0) ;
}

int ZoltanSpace::BaseLB::GetCommonParameter(char *key, double *d) 
{

  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;
  
  stringstream dummy ; dummy << props[ ZoltanName ] << ends ;
  dummy >> *d ;
  
  return(0) ;
}

/* Automigrate : hard-coded to false.
   Deterministic : false (default true)
   UseMachFile : hard-coded to false now
*/
int ZoltanSpace::BaseLB::SetCommonParameter(char *key, bool b)
{

  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;

  string charform; 
  charform = (b == true) ? "1" : "0" ;

  props[ ZoltanName ] = charform ;
  Zoltan_Set_Param( my_zz, const_cast<char *> (ZoltanName.c_str()), 
		    const_cast<char *> (charform.c_str()) );

  return(0) ;  
}

int ZoltanSpace::BaseLB::GetCommonParameter(char *key, bool *b)
{
  
  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;
 
  string ans = props[ ZoltanName ] ;
  *b = (ans == "1") ? true : false ;

 return(0) ;
}

/* NumGidEntries :  a number, default 1
   NumLidEntries :  a number, default 1
   ObjWtDim      :  a number {0 / 1}, default 0
   EdgeWtDim     :  a number {0 / 1}, default 0
   DebugLevel    :  a number > 0 ; default 1.
   DebugProc     :  a number > 0 ; default 0.
   CommWtDim     :  a number > 0 ; default 1
*/
int ZoltanSpace::BaseLB::SetCommonParameter(char *key, int i)
{

  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;

  stringstream dummy ; dummy << i << ends ;
  string charform; dummy >> charform ;

  props[ ZoltanName ] = charform ;
  Zoltan_Set_Param( my_zz, const_cast<char *> (ZoltanName.c_str()), 
		    const_cast<char *>((props[ ZoltanName ]).c_str()) );

  return(0) ;
}

int ZoltanSpace::BaseLB::GetCommonParameter(char *key, int *i)
{

  if (is_init == false) init() ;

  string inname(key) ;
  if ( translation.find( inname ) == translation.end() ) return(-1) ;

  string ZoltanName = translation[ inname ] ;
  
  stringstream dummy ; dummy << props[ ZoltanName ] << ends ;
  dummy >> *i ;

  return(0) ;
}

void ZoltanSpace::BaseLB::init()
{
  props["NUM_GID_ENTRIES"] = "1" ;
  props["NUM_LID_ENTRIES"] = "1" ;
  props["DEBUG_LEVEL"] = "1" ;
  props["DEBUG_PROCESSOR"] = "0" ;
  props["TIMER"] = "wall" ;
  props["OBJ_WEIGHT_DIM"] = "0" ;
  props["EDGE_WEIGHT_DIM"] = "0" ;
  props["DETERMINISTIC"] = "TRUE" ;
  props["USE_MACHINE_DESC"] = "0" ;
  props["MACHINE_DESC_FILE"] = "/etc/local/MachineDesc" ;
  props["DEBUG_MEMORY"] = "1";
  props["IMBALANCE_TOL"] = "1.2" ;
  props["COMM_WEIGHT_DIM"] = "1" ;

  translation["Timer"] = "Timer" ;
  translation["Deterministic"] = "DETERMINISTIC" ;
  translation["UseMachFile"] = "USE_MACHINE_DESC" ;
  translation["NumGidEntries"] = "NUM_GID_ENTRIES" ;
  translation["NumLidEntries"] = "NUM_LID_ENTRIES" ;
  translation["ObjWtDim"] = "OBJ_WEIGHT_DIM" ;
  translation["EdgeWtDim"] = "EDGE_WEIGHT_DIM" ;
  translation["DebugLevel"] = "DEBUG_LEVEL" ;
  translation["DebugProc"] = "DEBUG_PROCESSOR" ;
  translation["ImbalanceTolerance"] = "IMBALANCE_TOL" ;
  translation["CommWtDim"] = "COMM_WEIGHT_DIM" ;

  is_init = true ;
}

int ZoltanSpace::BaseLB::PrintCommonKeywordsToScreen()
{
  if (is_init == false ) init() ;

  cout << " These common-to-all-load-balancers keywords and values are : " << endl ;
  map<string, string>::iterator p ;
  for( p = translation.begin() ; p != translation.end(); p++ )
    cout << " Key = " << (*p).first << "   : value = " << props[(*p).second] << endl ;

  return(0) ;
}
