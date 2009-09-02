/*
  Implementation of the factory for load-balancers
  Jaideep Ray, SNL, 08/27/01
*/

#include <string>
#include "PartitionerFactory.h"

#include "RCB.h"
#include "RIB.h"
#include "HSFC.h"
#include "OctTree.h"
#include "ParMetis.h"

ZoltanSpace::PartitionerFactory_JR::PartitionerFactory_JR()
{
  is_init = false ;
  pEntityInfo = 0 ; pGeomInfo = 0 ; pEdgeInfo = 0 ; pTreeInfo = 0 ; pDataMig = 0 ;
}

ZoltanSpace::PartitionerFactory_JR::~PartitionerFactory_JR()
{
  vector< ::LoadPartitionerSpace::LoadBalancer *>::iterator p;
  for( p = my_lbs.begin(); p != my_lbs.end(); p++)
   {
     delete (*p) ; *p = 0 ;
   }
  
}

void ZoltanSpace::PartitionerFactory_JR::setServices( CONV_NS(Services) *svc)
{

  if ( svc == 0 ) 
  {
    if ( is_init == true )
    {
      pSvc->releasePort("EntityInfo") ;
      pSvc->releasePort("GeomInfo") ;
      pSvc->releasePort("EdgeInfo") ;
      pSvc->releasePort("TreeInfo") ;
      pSvc->releasePort("DataMig") ;
      pSvc->releasePort("cSvc") ;
    }
    return ;
  }
  pSvc = svc ;

  pSvc->addProvidesPort( this, pSvc->createPortInfo("LoadBalancerFactory",
						    "LoadBalancerFactoryPort",
						    0) ) ;
  pSvc->addProvidesPort( this, pSvc->createPortInfo("DataLoc",
					    "::LoadPartitionerSpace::DataLocationPort",
					    0) ) ;
  pSvc->registerUsesPort( svc->createPortInfo("EntityInfo", 
					      "::LoadPartitionerSpace::EntityInfo", 0) );
  pSvc->registerUsesPort( svc->createPortInfo("GeomInfo", 
					      "::LoadPartitionerSpace::GeomInfo", 0) );
  pSvc->registerUsesPort( svc->createPortInfo("EdgeInfo", 
					      "::LoadPartitionerSpace::EdgeInfo", 0) );
  pSvc->registerUsesPort( svc->createPortInfo("TreeInfo", 
					      "::LoadPartitionerSpace::TreeInfo", 0) );
  pSvc->registerUsesPort( svc->createPortInfo("DataMigrate", 
					      "::LoadPartitionerSpace::DataMigrationPort", 0) );
  
  // Configurable parameter port
  svc->registerUsesPort( svc->createPortInfo("cSvc", 
					     "gov.cca.ParameterPortFactoryService", 0) ) ;
  ConfigurableParameterFactory *cpf = 
    dynamic_cast<ConfigurableParameterFactory *>(svc->getPort("cSvc"));
  CHECKDC(cpf) ;
  pp = cpf->createConfigurableParameterPort();  CHECKDC(pp) ;
  svc->addProvidesPort(pp, svc->createPortInfo("CONFIG", "ConfigurableParameterPort", 0) );

  pIterate = new BoolParameter("IteratorAccess", 
			       "Is data acessed as lists or by iterating ?",
			       "Iterator access =", true) ;
  pp->setBatchTitle("Partitioner Factory Parameters");
  pp->addRequest( pIterate ) ;  
}

::LoadPartitionerSpace::LoadBalancer *
ZoltanSpace::PartitionerFactory_JR::CreateLB_VM_pair(MPI_Comm *A, char *name)
{
  if (is_init == false) { init() ; is_init = true ; }
  string aname(name) ;
  
  if (aname == "RCB")
  {
    ::ZoltanSpace::RCB_LB *a = new ::ZoltanSpace::RCB_LB(this, my_lbs.size(), A) ;
    my_lbs.push_back( a ) ;
    return( dynamic_cast< ::LoadPartitionerSpace::LoadBalancer *> (a) ) ;
  }
  else if (aname == "RIB")
  {
    ::ZoltanSpace::RIB_LB *a = new ::ZoltanSpace::RIB_LB(this, my_lbs.size(), A) ;
    my_lbs.push_back( a ) ;
    return( dynamic_cast< ::LoadPartitionerSpace::LoadBalancer *> (a) ) ;
  }
  else if (aname == "HSFC")
  {
    ::ZoltanSpace::HSFC_LB *a = new ::ZoltanSpace::HSFC_LB(this, my_lbs.size(), A) ;
    my_lbs.push_back( a ) ;
    return( dynamic_cast< ::LoadPartitionerSpace::LoadBalancer *> (a) ) ;
  }
  else if (aname == "OCTPART")
  {
    ::ZoltanSpace::OctTree_LB *a = new ::ZoltanSpace::OctTree_LB(this, 
								 my_lbs.size(), A) ;
    my_lbs.push_back( a ) ;
    return( dynamic_cast< ::LoadPartitionerSpace::LoadBalancer *> (a) ) ;
  }
  else if (aname == "PARMETIS")
  {
    ::ZoltanSpace::ParMetis_LB *a = new ::ZoltanSpace::ParMetis_LB(this, 
								   my_lbs.size(), A) ;
    my_lbs.push_back( a ) ;
    return( dynamic_cast< ::LoadPartitionerSpace::LoadBalancer *> (a) ) ;
  }
  else
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::CreateLB_VM_pair() "
	 << " load balancer of type " << aname << " not known. return with 0 " 
	 << endl ;
  }
  return(0) ;
}

void ZoltanSpace::PartitionerFactory_JR::Destroy_LB_VM_pair( 
                                           ::LoadPartitionerSpace:: LoadBalancer *p)
{
  if ( is_init == false ) return ;

  p->cleanup() ;  delete p ; 

  return ;
}

void ZoltanSpace::PartitionerFactory_JR::init()
{
  pEntityInfo = dynamic_cast< ::LoadPartitionerSpace::EntityInfo *> 
                 ( pSvc->getPort("EntityInfo") ) ;
  if (pEntityInfo == 0)
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::init()"
	 << " pEntityInfo is NULL. Caution !" << endl;
  }

  pGeomInfo = dynamic_cast< ::LoadPartitionerSpace::GeomInfo *> 
                 ( pSvc->getPort("GeomInfo") ) ;
  if (pGeomInfo == 0)
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::init()"
	 << " pGeomInfo is NULL. Caution !" << endl;
  }

  pEdgeInfo = dynamic_cast< ::LoadPartitionerSpace::EdgeInfo *> 
                 ( pSvc->getPort("EdgeInfo") ) ;
  if (pEdgeInfo == 0)
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::init()"
	 << " pEdgeInfo is NULL. Caution !" << endl;
  }

  pTreeInfo = dynamic_cast< ::LoadPartitionerSpace::TreeInfo *> 
                 ( pSvc->getPort("TreeInfo") ) ;
  if (pTreeInfo == 0)
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::init()"
	 << " pTreeInfo is NULL. Caution !" << endl;
  }

  pDataMig = dynamic_cast< ::LoadPartitionerSpace::DataMigrationPort *>
                  ( pSvc->getPort("DataMigrate") ) ;
  if ( pDataMig == 0)
  {
    cerr << " ZoltansSpace::PartitionerFactor_JR::init() "
	 << " pDataMig is NULL. Caution !" << endl ;
  }
  is_init = true ;
}
