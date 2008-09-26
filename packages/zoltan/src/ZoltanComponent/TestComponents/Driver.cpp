/*
  Implementation of the driver.
*/

#include "Driver.h"
#include "ports/MPIBorrow.h"

extern "C"
{
#include "dr_input_const.h"      /* where read_cmd_file, etc are */
#include "dr_err_const.h"        /* error report */

  //These are things that I do not use, but are needed by the library I
  // link to.
  int Debug_Driver = 1 ;
  int Number_Iterations = 1 ; 
  int Driver_Action = 1 ;
  int Debug_Chaco_Input = 0 ;
  int Chaco_In_Assign_Inv = 0 ;
  struct Test_Flags Test ;
  struct Output_Flags Output ;
}

ZoltanTestSpace::Driver::Driver()
{
  is_init = false ;
  pInpFile = new StringParameter("CmdFilename", " Command Filename = ", 
				 "Zoltan test input files", "zdrive.inp.rcb") ;

  ZoltanToMine[ "Timer" ] = "Timer" ;
  ZoltanToMine[ "DETERMINISTIC" ] = "Deterministic" ;
  ZoltanToMine[ "USE_MACHINE_DESC" ] = "UseMachFile" ;
  ZoltanToMine[ "NUM_GID_ENTRIES" ] = "NumGidEntries" ;
  ZoltanToMine[ "NUM_LID_ENTRIES" ] = "NumLidEntries" ;
  ZoltanToMine[ "OBJ_WEIGHT_DIM" ] = "ObjWtDim" ;
  ZoltanToMine[ "EDGE_WEIGHT_DIM" ] = "EdgeWtDim" ;
  ZoltanToMine[ "DEBUG_LEVEL" ] =  "DebugLevel" ;
  ZoltanToMine[ "DEBUG_PROCESSOR" ] = "DebugProc" ;
  ZoltanToMine[ "IMBALANCE_TOL" ] = "ImbalanceTolerance" ;
  ZoltanToMine[ "COMM_WEIGHT_DIM" ] = "CommWtDim" ;

  pio_info.dsk_list_cnt = -1;     pio_info.num_dsk_ctrlrs = -1 ;
  pio_info.pdsk_add_fact = -1;    pio_info.zeros = -1;
  pio_info.file_type = -1 ;       pio_info.init_dist_type = -1;
  pio_info.init_size = -1 ;       pio_info.init_dim = -1 ;
  pio_info.init_vwgt_dim = -1 ;
  pio_info.pdsk_root[0] = '\0';      pio_info.pdsk_subdir[0] = '\0' ;
  pio_info.pexo_fname[0] = '\0';  

  prob.method[0] = '\0';  prob.num_params = 0  ; prob.params = NULL ;

  // Initialize the global structs, Test and Output
  Test.DDirectory = 0 ;
  Test.Local_Parts = 0 ;
  Test.Drops = 0 ;
  Test.Multi_Callbacks = 0 ;
  Test.Gen_Files = 0 ;
  Test.Null_Lists = NONE ;

  Output.Gnuplot = 0 ;
  Output.Nemesis = 0 ;
  Output.Plot_Partition = 0 ;
  Output.Mesh_Info_File = 0 ;
}

ZoltanTestSpace::Driver::~Driver()
{
  delete pInpFile ;
  if (is_init == true)
  {
    delete (prob.params) ;
  }
}

void ZoltanTestSpace::Driver::setServices( CONV_NS(Services) *svc)
{

  if ( svc == NULL )
  {
    pSvc->releasePort("cSvc");  pSvc->releasePort("LBF") ; pSvc->releasePort("IO") ;
    pSvc->releasePort("mpi") ;  pSvc->releasePort("go") ; pSvc->releasePort("CONFIG") ;
    pSvc->releasePort("DataMig") ;

    return ;
  }

  this->pSvc = svc ;

  // Tell framework about the ports i use
  pSvc->registerUsesPort( svc->createPortInfo("LBF", 
				"::LoadPartitionerSpace::LoadBalancerFactory", 0) );
  pSvc->registerUsesPort( svc->createPortInfo("DataLoc", 
				   "::LoadPartitionerSpace::DataLocationPort", 0) ) ;
  pSvc->registerUsesPort( 
                       svc->createPortInfo("IO", "::ZoltanTestSpace::IOPort", 0) ) ;
  svc->registerUsesPort( svc->createPortInfo("mpi", "gov.cca.MPIBorrow",0) );

  // I present this port
  pSvc->addProvidesPort( this, svc->createPortInfo("go", "GoPort", 0) ) ;

  // Configurable parameter for cmd file
  svc->registerUsesPort( svc->createPortInfo("cSvc", 
				      "gov.cca.ParameterPortFactoryService", 0) ) ;
  ConfigurableParameterFactory *cpf = 
                 dynamic_cast<ConfigurableParameterFactory *>(svc->getPort("cSvc"));
  pp = cpf->createConfigurableParameterPort();   
  svc->addProvidesPort(pp, 
                    svc->createPortInfo("CONFIG", "ConfigurableParameterPort", 0) );

  pp->setBatchTitle("Driver parameters");
  pp->addRequest( pInpFile ) ;
}

void ZoltanTestSpace::Driver::init()
{
  int ans ;

  ans = read_cmd_file( pInpFile->value, &prob, &pio_info, NULL) ;  
  if ( ans != 1 ) 
  { 
    error_report("ZoltanTestSpace::Driver::init()", __LINE__, __FILE__) ;
    std::terminate() ;
  }

  ans = check_inp(&prob, &pio_info) ;  
  if ( ans != 1 ) 
  { 
    error_report("ZoltanTestSpace::Driver::init()", __LINE__, __FILE__) ;
    std::terminate() ;
  }

  CONV_NS(MPIBorrow) *pMPI_Port = dynamic_cast< CONV_NS(MPIBorrow) *> 
                                  (pSvc->getPort("mpi") );
  int key, tag;  comm = pMPI_Port->borrowComm(1, &tag, key) ;
  MPI_Comm_rank(comm, &Proc) ; MPI_Comm_size(comm, &Num_Procs) ;

  //JR_DEBUG IS this really needed ? since all the procs do read individually
  //  brdcst_cmd_info(Proc, &prob, &pio_info) ;

  /* now calculate where the file for this processor is */
  if(pio_info.dsk_list_cnt <= 0) 
  {
    if (pio_info.num_dsk_ctrlrs > 0) 
    {
      int ctrl_id = (Proc % pio_info.num_dsk_ctrlrs);
      pio_info.rdisk = ctrl_id + pio_info.pdsk_add_fact;
    }
  }
  else 
  {
    int ctrl_id = Proc % pio_info.dsk_list_cnt;
    pio_info.rdisk = pio_info.dsk_list[ctrl_id];
  }

  pLBF = dynamic_cast< ::LoadPartitionerSpace::LoadBalancerFactory * > 
    ( pSvc->getPort("LBF") ) ; CHECKDC(pLBF) ;

  pIO = dynamic_cast< ::ZoltanTestSpace::IOPort * > ( pSvc->getPort("IO") ) ;
  CHECKDC(pIO) ;
  pIO->SetCmdFilename( pInpFile->value ) ; // Tell the mesh where to read inputs from.

  pDLP = dynamic_cast< ::LoadPartitionerSpace::DataLocationPort * > 
                       ( pSvc->getPort("DataLoc") ) ; CHECKDC(pDLP) ;
  is_init = true ;
}

int ZoltanTestSpace::Driver::go()
{
  if (is_init == false) { init() ; is_init = true ; }

  // Create a load balancer
  std::string lb_name( prob.method ) ;
  std::transform( lb_name.begin(), lb_name.end(), lb_name.begin(), toupper) ;
  char *name = const_cast< char *> (lb_name.c_str()) ;
  ::LoadPartitionerSpace::LoadBalancer *pRCB = pLBF->CreateLB_VM_pair(&comm, name );

  // Configure load-balancer
  SetParams( pRCB ) ; 

  // Dump out grid
  pIO->ProduceGnuplotOutput( "Before" ) ;  
  pIO->ProduceFileOutput( "Before" ) ;  

  // Create a balanced partition
  bool did_it_change ;
  ::LoadPartitionerSpace::EntityList *pImport, *pExport ;
  int ans = pRCB->CreateBalancedPartition( &did_it_change, &pImport, &pExport) ;
  if ( ans != 0 ) 
  {
    std::cerr << " ZoltanTestSpace::Driver::go() load balance failed " << std::endl;
    return(-1) ;
  }
  else  // Write out stats about what needs to be moved
  {
    std::cout << "Proc = " << Proc << " did it change = " << did_it_change 
	 << " Import = " << pImport->GetListLength() 
	 << " Export = " << pExport->GetListLength() 
	 << std::endl ;
    MPI_Barrier(comm) ;
  }
  
  // Create a point and see which processors it gets assigned to.
  // This will work for RCB, RIB and Hilbert Space Filling Curves partitioners, others
  // return with -1.
  if (Proc == 0 ) 
  {
    std::cout << std::endl << std::endl << " ------ Point and box tests ------ " 
	      << std::endl << std::endl ;

    double coords[2] = { 0.0, 0.0} ;
    int goes_onto_this_proc ;
    ans = pRCB->IncrementalAssignPoint(coords, 2, &goes_onto_this_proc) ;
    if ( ans == 0 ) 
      std::cout << " The point [0,0] goes onto proc " << goes_onto_this_proc << std::endl ;
    
    // Create a box and see which procs it gets divided into
    // This will work for RCB, RIB and Hilbert Space Filling Curves partitioners, others
    // return with -1.
    double lbbc[2] = {0.0, 0.0} ;    double ubbc[2] = {1.0, 1.0} ;
    int nprocs, *proc_list ;
    ans = pRCB->IncrementalAssignBox(lbbc, ubbc, 2, &nprocs, &proc_list) ;
    if (ans == 0) 
    {
      std::cout << " The box {0,0} {1,1} goes onto procs : " ;
      for(int p = 0; p < nprocs; p++) std::cout << " " << proc_list[p] << " " ;
      std::cout << std::endl ;
    }
    std::cout  << std::endl << " ------------------------------------ " << std::endl ;
  }
  MPI_Barrier(comm) ;

  // Now move the entities to realize the balanced partition, then cleanup
  pImport->addRef() ; pExport->addRef() ;
  pDLP->MigrateEntities( pRCB, pExport, pImport) ;
  pDLP->Cleanup( pRCB ) ;
  
  // Dump out load-balanced grid
  pIO->ProduceGnuplotOutput( "After" ) ; 
  pIO->ProduceFileOutput( "After" ) ; 

  // I'm done with these
  pImport->deleteRef() ; pExport->deleteRef() ; 

  // remove any data or temporaries, prior to being used again
  pRCB->cleanup() ; 
  pIO->DeleteAllAndRestart() ;

  // Done with this load-balancer; kill it.
  pLBF->Destroy_LB_VM_pair( pRCB ) ;

  return(0) ;
}

int ZoltanTestSpace::Driver::SetParams( ::LoadPartitionerSpace::LoadBalancer *pLB)
{
  for(int i = 0;  i < prob.num_params; i++)
  {
    std::string zoltan_name(prob.params[i][0]) ;
    transform( zoltan_name.begin(), zoltan_name.end(), zoltan_name.begin(), toupper);

    std::string my_name;
    if ( ZoltanToMine.find(zoltan_name) != ZoltanToMine.end() )
      my_name = ZoltanToMine[ zoltan_name ] ;
    else
      my_name = zoltan_name ;
    char *cc = const_cast< char * > (my_name.c_str()) ;

    if ( (ZoltanToMine.find( zoltan_name  ) == ZoltanToMine.end() ) &&
	 (zoltan_name != "PARMETIS_METHOD") )
    {
      std::cerr << " ZoltanTestSpace::Driver::SetParams() Did not find "
	   << " a parameter named " << zoltan_name << "; assummed integer value "
	   << std::endl ;
      /*
	OK, take a risk. If I haven't determined the translation, chances
	are that the keyword coming in from the input file is the same as
	accepted by the load-balancer. Further, the value corresponding to
	the keyword is usually an integer, so transfer it as such, by default.
	but warn before hand.
      */
    }

    if ( my_name == "Timer" || my_name == "PARMETIS_METHOD") 
      pLB->SetParameter( cc, prob.params[i][1] ) ;
    else if ( my_name == "ImbalanceTolerance" ) 
    {
      std::stringstream G ; G << prob.params[i][1] << std::ends;
      double dd ; G >> dd ;
      char *cc = const_cast< char * > (my_name.c_str()) ;
      pLB->SetParameter( cc, dd) ;
    }
    else if ( my_name == "Deterministic" || my_name == "UseMachFile" ) 
    {
      bool what = false ;
      if ( strcmp(prob.params[i][1], "1") == 0 ) what = true ;
      char *cc = const_cast< char * > ( my_name.c_str() ) ;
      pLB->SetParameter(cc, what) ;
    }
    else
    {
      // By default assume all values are ints
      std::stringstream G ; G << prob.params[i][1] << std::ends;
      int ii ; G >> ii ; char *cc = const_cast< char * > (my_name.c_str()) ;
      int ans = pLB->SetParameter(cc, ii) ; 
    }

  } // Loop over all problem parameters
  return(0) ;
}
      
      
    
