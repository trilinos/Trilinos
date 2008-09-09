/*
  Implementation of the Zoltan-free part of Mesh
*/

#include "Mesh.h"
#include "mpi.h"
#include "ports/MPIBorrow.h"

extern "C"
{
#include "dr_elem_util_const.h"  /* this is where free_mesh_arrays is defined */
#include "dr_output_const.h" /* this is where output_gnu  is */
#include "dr_err_const.h"

  // Some globals that are need for the test suite in ../../Zoltan/driver
  int Debug_Driver = 1 ;
  int Number_Iterations = 1 ; 
  int Driver_Action = 1 ;
  int Debug_Chaco_Input = 0 ;
  int Chaco_In_Assign_Inv = 0 ;
  struct Test_Flags Test ;
  struct Output_Flags Output ;

  double Total_Partition_Time = 0.0 ;
}

ZoltanTestSpace::Mesh::Mesh()
{
  is_init = false ;
  cmd_filename_set = false ;
  initialize_mesh() ;

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

ZoltanTestSpace::Mesh::~Mesh()
{
  if (is_init == true)
  {
    free_mesh_arrays( &mesh ) ;
    delete (prob.params) ;
    pSvc->unregisterUsesPort("cSvc");
  }
}

void ZoltanTestSpace::Mesh::setServices( CONV_NS(Services) *svc)
{
  this->pSvc = svc ;

  if ( svc == NULL )
  {
    pSvc->unregisterUsesPort("cSvc");
    return ;
  }

  // Tell framework about the ports i provide
  pSvc->addProvidesPort(this, pSvc->createPortInfo("EntityInfo", 
			      "::LoadBalancerSpace::EntityInfoPort", 0) ) ;
  pSvc->addProvidesPort(this, pSvc->createPortInfo("GeomInfo", 
			      "::LoadBalancerSpace::GeomInfoPort", 0) ) ;
  pSvc->addProvidesPort(this, pSvc->createPortInfo("EdgeInfo", 
			      "::LoadBalancerSpace::EdgeInfoPort", 0) ) ;
  pSvc->addProvidesPort(this, pSvc->createPortInfo("IO", 
			      "::ZoltanTestSpace::IOPort", 0) ) ;
  pSvc->addProvidesPort(this, pSvc->createPortInfo("DataMig", 
			      "::LoadBalancerSpace::DataMigrationPort", 0) ) ;

  // Notify that i use the mpi port to get a communicator.
  svc->registerUsesPort( svc->createPortInfo("mpi", "gov.cca.MPIBorrow",0) );
}

void ZoltanTestSpace::Mesh::initialize_mesh()
{
  /* Initializes mesh variables */

  mesh.data_type = MESH ;
  mesh.num_elems = mesh.num_nodes
                  = mesh.num_dims
                  = mesh.num_el_blks
                  = mesh.num_node_sets
                  = mesh.num_side_sets
                  = mesh.necmap
                  = mesh.elem_array_len
                  = mesh.nhedges
                  = mesh.vwgt_dim
                  = mesh.ewgt_dim
                  = mesh.hewgt_dim
                  = 0;
  mesh.eb_names                 = NULL;
  mesh.eb_ids                   = NULL;
  mesh.eb_cnts                  = NULL;
  mesh.eb_nnodes                = NULL;
  mesh.eb_nattrs                = NULL;
  mesh.ecmap_id                 = NULL;
  mesh.ecmap_cnt                = NULL;
  mesh.ecmap_elemids            = NULL;
  mesh.ecmap_sideids            = NULL;
  mesh.ecmap_neighids           = NULL;
  mesh.elements                 = NULL;
  mesh.hindex                   = mesh.hvertex = NULL ;
  mesh.hewgts                   = NULL ;
}

void ZoltanTestSpace::Mesh::init()
{
  if (cmd_filename_set == false)
  {
    std::cerr << "You need to set the command filename via the IOPort. Exiting"
	 << std::endl ;
    exit(-1) ;
  }

  int ans ;
  char *cc = const_cast< char * >  ( cmd_filename.c_str() ) ;
  ans = read_cmd_file( cc, &prob, &pio_info, NULL) ;
  if ( ans != 1 ) 
  { 
    std::cerr << " ZoltanTestSpace::Mesh::init() Error report 1 " << std::endl ;
    error_report(Proc); 
    std::cerr << "Exit at line "<< __LINE__ << " in file " << __FILE__ << std::endl ;
    std::terminate() ;

  }

  ans = check_inp(&prob, &pio_info) ;
  if ( ans != 1 ) 
  { 
    std::cerr << " ZoltanTestSpace::Mesh::init() Error report 2 " << std::endl ;
    error_report(Proc); 
    std::cerr << "Exit at line "<< __LINE__ << " in file " << __FILE__ << std::endl ;
    std::terminate() ;
  }

  CONV_NS(MPIBorrow) *pMPI_Port = dynamic_cast< CONV_NS(MPIBorrow) *> 
                                              (pSvc->getPort("mpi") );
  int key, tag;
  MPI_Comm comm = pMPI_Port->borrowComm(1, &tag, key) ;
  MPI_Comm_rank(comm, &Proc) ;
  MPI_Comm_size(comm, &Num_Procs) ;

  brdcst_cmd_info(Proc, &prob, &pio_info, &mesh) ;

  read_mesh() ;
  
  is_init = true ;
}

void ZoltanTestSpace::Mesh::read_mesh()
{
  if (pio_info.file_type == CHACO_FILE) 
  {
    if (!read_chaco_file(Proc, Num_Procs, &prob, &pio_info, &mesh)) 
    {
      std::cerr << "  ZoltanTestSpace::Mesh::read_mesh() " 
	   << "fatal: Error returned from read_chaco_mesh" << std::endl ;
    }
  }
  else if (pio_info.file_type == NEMESIS_FILE) 
  {
    if (!read_exoII_file(Proc, Num_Procs, &prob, &pio_info, &mesh)) 
    {
      std::cerr << " ZoltanTestSpace::Mesh::read_mesh() "
	   << " fatal: Error returned from read_exoII_mesh" << std::endl;
    }
  }
  else if (pio_info.file_type == HYPERGRAPH_FILE) {
    if (!read_hypergraph_file(Proc, Num_Procs, &prob, &pio_info, &mesh)) 
    {
      std::cerr << " ZoltanTestSpace::Mesh::read_mesh() "
		<< " fatal: Error returned from read_hypergraph_file" << std::endl;
    }
  }
  else if (pio_info.file_type == NO_FILE) {
    if (!create_random_input(Proc, Num_Procs, &prob, &pio_info, &mesh)) 
    {
      std::cerr << " ZoltanTestSpace::Mesh::read_mesh() "
		<< " fatal: Error returned from create_random_input" << std::endl ;
    }
  }
  else 
  {
    std::cerr << " ZoltanTestSpace::Mesh::read_mesh() "
	      << " fatal: Invalid file type" << std::endl ;
  }

}
 
int ZoltanTestSpace::Mesh::DeleteAllAndRestart()
{
  if (is_init == false) return(0) ;

  free_mesh_arrays( &mesh ) ;
  initialize_mesh() ;

  is_init = false ;
  return(0) ;
}

int ZoltanTestSpace::Mesh::ProduceGnuplotOutput(char *tag)
{
  if (is_init == false) { init() ; is_init = true ; } 

  char *cc = const_cast< char * > ( cmd_filename.c_str() ) ;
  output_gnu(cc, tag, Proc, Num_Procs, &prob, &pio_info, &mesh);
}

int ZoltanTestSpace::Mesh::SetCmdFilename(char *c)
{
  cmd_filename.append( c ) ;
  cmd_filename_set = true ;

  return(0) ;
}

int ZoltanTestSpace::Mesh::ProduceFileOutput(char *tag)
{
  if (is_init == false) { init() ; is_init = true ; } 

  char *cc = const_cast< char * > ( cmd_filename.c_str() ) ;
  output_results(cc, tag, Proc, Num_Procs, &prob, &pio_info, &mesh);
}
