/*
  The header for the Mesh component. This will have the ability to read
  Chaco & Nemesis mesh files in ../ZoltanTests. It will make extensive use
  of the code in ../../Zoltan/driver and ../../Zoltan/ch
*/

#ifndef DriverHSeen
#define DriverHSeen

// CCA standard definition of Port and Component
#include "cca.h"
#include "stdPorts.h"

//Sandia's own definition of the (user-) configurable parameters port
#include "jc++/jc++.h"
#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "port/supportInterfaces.h"
#include "dc/port/DefaultParameterPort.h"

//Ports I use
#include "IOPort.h"
#include "EntityList.h"
#include "LoadBalancer.h"
#include "LoadBalancerFactoryPort.h"
#include "DataLocationPort.h"

// C++ 
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>

// Specific to the driver that Zoltan has in ../../Zoltan/zdrive/
// Contains the definition of MESH_INFO, PARIO_INFO, PROB_INFO
#include "dr_const.h"
#include "dr_input_const.h"

#include "mpi.h"

// namespace conversion
#include "CONV_NS.h"

namespace ZoltanTestSpace
{
  class Driver : public virtual CONV_NS(Component), 
                 public virtual CONV_NS(GoPort)
  {
    public :

    Driver() ;

    virtual ~Driver() ;

    virtual void setServices( CONV_NS(Services) *svc) ;

    virtual int go() ;

    private :

    void init() ;
    
    int SetParams( ::LoadPartitionerSpace::LoadBalancer *pLB ) ;

    void error_report( char *mesg, int line, char *file) 
    {
      std::cerr << mesg << " on line " << line << " from " << file << std::endl ;
    }

    PARIO_INFO pio_info ;
    PROB_INFO prob ;

    bool is_init ;
    int Proc, Num_Procs, num_gid_entries, num_lid_entries ;
    StringParameter *pInpFile ;
    CONV_NS(Services) *pSvc ;
    
    // Ports i use 
    ::LoadPartitionerSpace::LoadBalancerFactory *pLBF;
    ::ZoltanTestSpace::IOPort *pIO;
    ::LoadPartitionerSpace::DataLocationPort *pDLP ;

    std::map< std::string, std::string> ZoltanToMine ;
    ConfigurableParameterPort *pp ;
    MPI_Comm comm ;
  } ;
} ;
#endif
      
    
