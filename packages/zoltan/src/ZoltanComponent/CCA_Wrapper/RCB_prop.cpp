#include "RCB.h"
#include <sstream>

int ZoltanSpace::RCB_LB::SetParameter(char *key, char *val)
{

  if (is_init == false) init() ;
  
  string name(key) ; int err = -1 ;
  if ( prop.find(name) == prop.end() ) err = BaseLB::SetCommonParameter(key, val) ;

  if ( err != 0)
  {
    string value(val) ;
    prop[ name ] = value ;
    Zoltan_Set_Param(BaseLB::my_zz, const_cast<char *> (name.c_str()), val) ;
  }
  return(0) ;
}

int ZoltanSpace::RCB_LB::GetParameter(char *key, char **val)
{

  if (is_init == false) init() ;

  string name(key) ;
  if ( prop.find(name) == prop.end() ) 
    return( BaseLB::GetCommonParameter(key, val) ) ;

  *val = const_cast< char * > ( ( prop[ name ] ).c_str() ) ;

  return(0) ;
}

int ZoltanSpace::RCB_LB::SetParameter(char *key, double d)
{

  if (is_init == false) init() ;

  string name(key) ; int err = -1 ;
  if ( prop.find(name) == prop.end() ) err = BaseLB::SetCommonParameter(key, d) ;
  
  if ( err != 0 )
  {
    stringstream dummy ;
    dummy << d << ends ;
    string G; dummy >> G ;
    
    prop[name] = G ;
    Zoltan_Set_Param(BaseLB::my_zz, 
		     const_cast< char *> (name.c_str()), const_cast<char *> (G.c_str()) );
  }
  return(0);
}
 
int ZoltanSpace::RCB_LB::GetParameter(char *key, double *d)
{

  if (is_init == false) init() ;

  string name(key) ;
  if ( prop.find(name) == prop.end() ) 
    return( BaseLB::GetCommonParameter(key, d) ) ;

  stringstream dummy ;
  dummy << prop[ name ] << ends ;
  dummy >> *d ;

  return(0);
}
    
int ZoltanSpace::RCB_LB::SetParameter(char *key, bool b)
{

  if (is_init == false) init() ;

  string name(key) ; int err = -1 ;
  if ( prop.find(name) == prop.end() ) err = BaseLB::SetCommonParameter(key, b) ;

  if ( err != 0 ) 
  {
    string value = (b == true) ? "1" : "0" ;
    prop[ name ] = value ;
    Zoltan_Set_Param(BaseLB::my_zz, key, const_cast<char *>(value.c_str()));
  }
  return(0);
}

int ZoltanSpace::RCB_LB::GetParameter(char *key, bool *b)
{

  if (is_init == false) init() ;

  string name(key) ;
  if ( prop.find(name) == prop.end() )
    return( BaseLB::GetCommonParameter(key, b) ) ;

  string G = prop[ name ] ;
  *b = ( G == "1" ) ? true : false ;

  return(0);
}
    
int ZoltanSpace::RCB_LB::SetParameter(char *key, int i)
{

  if (is_init == false) init() ;

  string name(key) ; int err = -1 ;
  if ( prop.find(name) == prop.end() ) err = BaseLB::SetCommonParameter(key, i) ;

  if ( err != 0 )
  {
    stringstream dummy ;
    dummy << i << ends;
    string g ;  dummy >> g ;
    
    prop[name] = g ;
    Zoltan_Set_Param(BaseLB::my_zz, const_cast<char *>(name.c_str()), 
		     const_cast<char *>(g.c_str()) );
  }
  return(0);
}

int ZoltanSpace::RCB_LB::GetParameter(char *key, int *i)
{

  if (is_init == false) init() ;

  string name(key) ;
  if ( prop.find(name) == prop.end() )
    return( BaseLB::GetCommonParameter(key, i) ) ;

  string G = prop[name] ;
  stringstream dummy ; dummy << G << ends;
  dummy >> *i ;

  return(0);
}

int ZoltanSpace::RCB_LB::PrintKeywordsToScreen()
{

  if (is_init == false) init() ;

  if (rank == master)
  {
    // Specific to me.
    cout << "The property keywords and values specific to this " 
	 << myname << " load-balancer are :" << endl ;

    ::map<string, string>::iterator p;
    for(p = prop.begin() ; p != prop.end(); p++)
      cout << " Key = " << p->first << "  : value " << p->second << endl ;

    // Print the common keywords and value
    BaseLB::PrintCommonKeywordsToScreen() ;
  }
  return(0) ;
}
