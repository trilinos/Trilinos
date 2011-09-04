
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include <stk_rebalance/ZoltanPartition.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <zoltan.h>

#include <Teuchos_ParameterList.hpp>

using namespace std;
using namespace stk;
using namespace stk::rebalance;

#define STK_GEOMDECOMP_DEBUG 0

namespace {


double static_zoltan_version(const double v=0) {
  static double version=0;
  if (v) { version=v;}
  return version;
}

inline unsigned num_gid_entries() {
  static const unsigned n=2;
  return n;
}

inline unsigned num_lid_entries() {
  static const unsigned n=2;
  return n;
}
inline unsigned wdim() {
  static const unsigned n=1;
  return n;
}

inline void convert_param_to_string(const Parameters &from,
                                    vector < pair<std::string, std::string> > &to)
{
  Parameters::ConstIterator
    from_iter  = from.begin(),
    from_end   = from.end();

  for (; from_iter != from_end; ++from_iter) {
    std::string s;
    const std::string name  = from.name(from_iter);
    const std::string entry = from.entry(from_iter).getValue(&s);
    std::pair<std::string,std::string> p(name,entry);

    to.push_back(p);
  }
}

inline void fill_parameters (const char *T[][2],
                             const int  i,
                             Parameters &Entry)
{
  for (int j=0; j<i; ++j) Entry.set(T[j][0], T[j][1]);
}

void fill_name_conversion( Parameters & name_conversion )
{

  Parameters & general = name_conversion.sublist("General");

  general.set("LOAD BALANCING METHOD"      , "LB_METHOD");
  general.set("ZOLTAN DEBUG LEVEL"         , "DEBUG_LEVEL");
  general.set("DEBUG PROCESSOR NUMBER"     , "DEBUG_PROCESSOR");
  general.set("TIMER"                      , "TIMER");
  general.set("DETERMINISTIC DECOMPOSITION", "DETERMINISTIC");
  general.set("DEBUG MEMORY"               , "DEBUG_MEMORY");
  general.set("IMBALANCE TOLERANCE"        , "IMBALANCE_TOL");
  general.set("RENUMBER PARTITIONS"        , "REMAP");
  general.set("KEEP CUTS"                  , "KEEP_CUTS");
  general.set("REUSE CUTS"                 , "RCB_REUSE");
  general.set("RCB RECOMPUTE BOX"          , "RCB_RECOMPUTE_BOX");
  general.set("CHECK GEOMETRY"             , "CHECK_GEOM");
  general.set("LOCK RCB DIRECTIONS"        , "RCB_LOCK_DIRECTIONS");
  general.set("SET RCB DIRECTIONS"         , "RCB_SET_DIRECTIONS");
  general.set("RCB MAX ASPECT RATIO"       , "RCB_MAX_ASPECT_RATIO");
  general.set("RECTILINEAR RCB BLOCKS"     , "RCB_RECTILINEAR_BLOCKS");
  general.set("OCTREE DIMENSION"           , "OCT_DIM");
  general.set("OCTREE METHOD"              , "OCT_METHOD");
  general.set("OCTREE MIN ENTITIES"         , "OCT_MINOBJECTS");
  general.set("OCTREE MAX ENTITIES"         , "OCT_MAXOBJECTS");
//  // These values are never changed, but must
  // be set so default values work correctly.
  general.set("NUMBER GLOBAL ID ENTRIES"   , "NUM_GID_ENTRIES");
  general.set("NUMBER LOCAL ID ENTRIES"    , "NUM_LID_ENTRIES");
  general.set("ENTITY WEIGHT DIMENSION"    , "OBJ_WEIGHT_DIM");
  general.set("RETURN LISTS"               , "RETURN_LISTS");
  general.set("AUTOMATIC MIGRATION"        , "AUTO_MIGRATE");
  general.set("DISTANCE"                   , "DISTANCE");

  Parameters & rcb = name_conversion.sublist("0");
  rcb.set("OVER ALLOCATE MEMORY"       , "RCB_OVERALLOC");
  rcb.set("ALGORITHM DEBUG LEVEL"      , "RCB_OUTPUT_LEVEL");

  Parameters & rib = name_conversion.sublist("1");
  rib.set("OVER ALLOCATE MEMORY"       , "RIB_OVERALLOC");
  rib.set("ALGORITHM DEBUG LEVEL"      , "RIB_OUTPUT_LEVEL");

  Parameters & hsfc = name_conversion.sublist("2");
  hsfc.set("OVER ALLOCATE MEMORY"       , "");
  hsfc.set("ALGORITHM DEBUG LEVEL"      , "" );

  Parameters & oct = name_conversion.sublist("3");
  oct.set("OVER ALLOCATE MEMORY"       , "");
  oct.set("ALGORITHM DEBUG LEVEL"      , "OCT_OUTPUT_LEVEL");

  Parameters & graph = name_conversion.sublist("4");
  graph.set("OVER ALLOCATE MEMORY"       , "");
  graph.set("ALGORITHM DEBUG LEVEL"      , "" );
}


void fill_value_conversion( Parameters & value_conversion )
{
  Parameters & lb_method = value_conversion.sublist("LOAD BALANCING METHOD");
  lb_method.set("0"   , "RCB");
  lb_method.set("1"   , "RIB");
  lb_method.set("2"   , "HSFC");
  lb_method.set("3"   , "OCTPART");
  lb_method.set("4"   , "GRAPH");

  Parameters & timer = value_conversion.sublist("TIMER");
  timer.set("0"   , "WALL");
  timer.set("1"   , "CPU");

}

void fill_default_values( Parameters & values )
{
  Parameters & default_values = values; //values.sublist("General");

  default_values.set("LOAD BALANCING METHOD"      , "0");
  default_values.set("RENUMBER PARTITIONS"        , "1");
  default_values.set("ZOLTAN DEBUG LEVEL"         , "0");
  default_values.set("TIMER"                      , "0");
  default_values.set("DETERMINISTIC DECOMPOSITION", "1");
  default_values.set("DEBUG MEMORY"               , "1");
  default_values.set("IMBALANCE TOLERANCE"        , "1.1");
  default_values.set("KEEP CUTS"                  , "1");
  default_values.set("REUSE CUTS"                 , "1");
  default_values.set("OVER ALLOCATE MEMORY"       , "1.1");
  default_values.set("ALGORITHM DEBUG LEVEL"      , "0");
//  default_values.set("OCTREE MIN ENTITIES"        , "1");
//  default_values.set("OCTREE MAX ENTITIES"        , "1");
  default_values.set("NUMBER GLOBAL ID ENTRIES"   , "2");
  default_values.set("NUMBER LOCAL ID ENTRIES"    , "2");
  default_values.set("ENTITY WEIGHT DIMENSION"    , "1");
  default_values.set("RETURN LISTS"               , "EXPORT");
}



#if STK_GEOMDECOMP_DEBUG>=2
void debug_print_decomp_export(Zoltan   *zoltan,
                               Zoltan *zoltan_id )
{
  int i;
  int           num_export   = zoltan->Num_Exported();
  ZOLTAN_ID_PTR export_lids  = zoltan->Export_Local_IDs();
  ZOLTAN_ID_PTR export_gids  = zoltan->Export_Global_IDs();
  int*          export_procs = zoltan->Export_Proc_IDs();
  int           Z_LID_SIZE   = zoltan->Num_Lid_Entries();
  int           Z_GID_SIZE   = zoltan->Num_Gid_Entries();

  if ( export_gids!=NULL && export_lids!=NULL && export_procs!=NULL ) {
    Env::output() << ": Zoltan RCB EXPORTS" << std::endl;
    for ( i = 0; i < num_export; i++ ) {
      Env::output()
        << "  " << i
        << ":  GID = "
        << "T" << zoltan_id->Type( &export_gids[i*Z_GID_SIZE]) << "  "
        << "I" << zoltan_id->Index(&export_gids[i*Z_GID_SIZE]) << "  "
        << "P" << zoltan_id->Proc( &export_gids[i*Z_GID_SIZE]) << "  "
        << "    LID = "
        << "T" << zoltan_id->Type( &export_lids[i*Z_LID_SIZE]) << "  "
        << "I" << zoltan_id->Index(&export_lids[i*Z_LID_SIZE]) << "  "
        << "  EXPORT_PROC_ID = "
        << export_procs[i]
        << std::endl;
    }
    Env::output_flush();
  }
}
#endif



extern "C" {
  int Callback_Num_Elements( void *data, int *ierr );
  void Callback_Element_List( void *data,
                             int Num_gid_entries,
                             int Num_lid_entries,
                             ZOLTAN_ID_PTR global_ids,
                             ZOLTAN_ID_PTR local_ids,     //{Iterator or index}
                             int wdim,
                             float *weights,
                             int *ierr );
  int Callback_Num_Dimensions( void *data, int *ierr );
  void Callback_Centroid_Coord( void *data,
                                int Num_gid_entries,
                                int Num_lid_entries,
                                ZOLTAN_ID_PTR global_id,
                                ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                                double *geom,
                                int *ierr );
  int Callback_Num_Edges( void *data,
                          int Num_gid_entries,
                          int Num_lid_entries,
                          ZOLTAN_ID_PTR global_id,
                          ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                          int *ierr );
  void Callback_Edge_List( void *data,
                           int Num_gid_entries,
                           int Num_lid_entries,
                           ZOLTAN_ID_PTR global_id,
                           ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                           ZOLTAN_ID_PTR nbor_global_id,
                           int *nbor_procs,
                           int wgt_dim,
                           float *ewgts,
                           int *ierr );
}


int Callback_Num_Elements( void *data, int *ierr )
{
  if ( data == NULL ) {
    *ierr = ZOLTAN_FATAL;  // Set FATAL Zoltan error flag
    return 0;
  }
  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan       *zdata = dynamic_cast<stk::rebalance::Zoltan*>(gdata);

  if (zdata == NULL ) {
    *ierr = ZOLTAN_FATAL;  // Set FATAL Zoltan error flag
    return 0;
  }

  int ne = zdata->num_elems();

  *ierr = ZOLTAN_OK;
  return ne;
}

void Callback_Element_List( void *data,
                            int Num_gid_entries,
                            int Num_lid_entries,
                            ZOLTAN_ID_PTR global_ids,
                            ZOLTAN_ID_PTR local_ids,     //{Iterator or index}
                            int weightdim,
                            float *weights,
                            int *ierr )
{
  if (!data) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return;
  }

  stk::rebalance::GeomDecomp   *gdata =  static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan       *zdata = dynamic_cast<stk::rebalance::Zoltan*>     (gdata);

  if (!zdata) {
    *ierr = ZOLTAN_FATAL;           // Set FATAL Zoltan error flag
    return;
  }

  unsigned k=0, l=0;
  const unsigned num_local_ids = zdata->num_moid();
  for (unsigned j=0; j<num_local_ids; ++j) {
    local_ids [k] = j;
    global_ids[k] = 0; // region_id
    ++k;
    local_ids [k] = 0;
    global_ids[k] = static_cast<ZOLTAN_ID_TYPE>(zdata->globalID(j));
    ++k;
    if (weightdim) {
      weights[l++] = zdata->entity_weight(j);
    } else {
      ++l;
    }
  }
  *ierr = ZOLTAN_OK;
  return;

}

int Callback_Num_Dimensions( void *data, int *ierr )
{
  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  const int  nd = zdata->spatial_dimension();

  *ierr = ZOLTAN_OK;
  return nd;

}

void Callback_Centroid_Coord( void *data,
                                     int Num_gid_entries,
                                     int Num_lid_entries,
                                     ZOLTAN_ID_PTR global_id,
                                     ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                                     double *geom,
                                     int *ierr )
{
  std::vector<double> temp(3,0.0);

  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }

  int lid = local_id[  0 ]; // Local Element ID

  const mesh::Entity & target_entity = * zdata->mesh_entity( lid );
  const GeomDecomp::VectorField & coor = * zdata->entity_coord_ref();
  const unsigned                        nd   =  zdata->spatial_dimension();

  /*
   * Obtain the centroid coordinates of the element by averaging all
   * the nodal coordinates of the nodes associated with the element.
   * Use GeomDecomp friend function, ent( , , )
   *
   */

  rebalance::GeomDecomp::entity_to_point( target_entity, coor, temp );

  for (size_t i=0 ; i < nd ; i++ ) geom[ i ] = (double) temp[ i ];

  *ierr = ZOLTAN_OK;
}



void getNeighbors( const mesh::Entity & entity,
                   std::set<const mesh::Entity*> & nodes )
{
  stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(entity);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  nodes.clear();

  mesh::PairIterRelation iElem = entity.relations(element_rank);

  for ( ; iElem.first != iElem.second; ++iElem.first ) {
    mesh::Entity * elem = iElem.first->entity();
    mesh::PairIterRelation iNode = elem->relations(fem_meta.node_rank());
    for ( ; iNode.first != iNode.second; ++iNode.first ) {
      mesh::Entity * node = iNode.first->entity();
      if (&entity != node) nodes.insert( node );
    }
  }
}

int numEdges( const mesh::Entity & entity ) {

  std::set<const mesh::Entity*> nodes;

  getNeighbors( entity, nodes );
  return nodes.size();
}

int Callback_Num_Edges( void *data,
                        int Num_gid_entries,
                        int Num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,  //{Iterator or index}
                        int *ierr )
{
  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  *ierr = ZOLTAN_OK;

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  int lid = local_id[  0 ]; // Local Element ID

  const mesh::Entity & target_entity = * zdata->mesh_entity( lid );

  return  numEdges( target_entity );

}

void Callback_Edge_List( void *data,
                         int Num_gid_entries,
                         int Num_lid_entries,
                         ZOLTAN_ID_PTR global_id,
                         ZOLTAN_ID_PTR local_id,  //{Iterator or index}

                         // Returned
                         ZOLTAN_ID_PTR nbor_global_id, //owning processor
                         int *nbor_procs,
                         int wgt_dim,
                         float *ewgts,

                         int *ierr )
{
  // (from include Fmwk_Sierra_Zoltan_Defines.h:)
  // ( Num_gid_entries = ZOLTAN_GID_SIZE 2 )
  // ( Num_lid_entries = ZOLTAN_LID_SIZE 2 )

  *ierr = ZOLTAN_OK;

  if ( !data ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }
  stk::rebalance::GeomDecomp *gdata = static_cast<stk::rebalance::GeomDecomp*>  (data);
  stk::rebalance::Zoltan     *zdata = dynamic_cast<stk::rebalance::Zoltan*>    (gdata);

  if ( !zdata ) {
    *ierr = ZOLTAN_FATAL;
    return ;
  }

  int lid = local_id[  0 ]; // Local Node ID

  const mesh::Entity & target_entity = * zdata->mesh_entity( lid );

  std::set<const mesh::Entity*> nodes;
  getNeighbors( target_entity, nodes );

  int counter(0);
  for ( std::set<const mesh::Entity*>::iterator i = nodes.begin();
        i != nodes.end(); ++i ) {
    nbor_global_id[counter*2+0] = 0; // region_id

    const mesh::Entity & n = **i;
    nbor_global_id[counter*2+1] = n.key().id();

    nbor_procs[counter] = n.owner_rank();

    if ( wgt_dim ) {
      ewgts[counter] = 1;
    }
    ++counter;
  }
}

}



const std::string Zoltan::m_zoltanparametersname_="Zoltan_Parameters";
const std::string Zoltan::m_defaultparametersname_="DEFAULT";


const std::string Zoltan::zoltan_parameters_name()
{
  return m_zoltanparametersname_;
}

const std::string Zoltan::default_parameters_name()
{
  return m_defaultparametersname_;
}

void Zoltan::init_default_parameters()
{
  fill_default_values(m_default_parameters_);
}


static Parameters *Name_Conversion =NULL;
static Parameters *Value_Conversion=NULL;

double Zoltan::zoltan_version()   const { return static_zoltan_version();  }

//: ===========
//: Constructor
//: ===========

namespace {
void merge_parameters(std::vector <std::pair<std::string, std::string> > &str_zoltan_params,
                      const Parameters &Zoltan_Params) {
  Parameters Merged_Zoltan_Params   ;
  Parameters Converted_Zoltan_Params;

  rebalance::Zoltan::merge_default_values (Zoltan_Params,
                                      Merged_Zoltan_Params);

  rebalance::Zoltan::convert_names_and_values(Merged_Zoltan_Params,
                                         Converted_Zoltan_Params);

  convert_param_to_string (Converted_Zoltan_Params,
                           str_zoltan_params);
  return;
}
}

Zoltan::Zoltan(ParallelMachine pm, const unsigned ndim, Parameters & rebal_region_parameters, const std::string parameters_name) :
  GeomDecomp(pm),
  m_zoltan_id_(NULL),
  m_spatial_dimension_(ndim),
  m_total_number_entities_(0)
{
  /* Determine if the default set of parameters already exists. */
  if( !rebal_region_parameters.isSublist(default_parameters_name()) )
  {
    init_default_parameters();
    rebal_region_parameters.sublist(default_parameters_name()) = m_default_parameters_;
  }

  /* If name is empty, use default values */
  std::string default_name =
    (parameters_name.empty()) ? default_parameters_name() : parameters_name ;

  if( !rebal_region_parameters.isSublist(default_name) ) {
    throw std::runtime_error("The Zoltan parameter set '" + default_name + "' does not exist.");
  }
  const Parameters & zoltan_params = rebal_region_parameters.sublist(default_name);

  std::vector <std::pair<std::string, std::string> > str_zoltan_params;
  merge_parameters(str_zoltan_params, zoltan_params);

  init(str_zoltan_params);
}

void
Zoltan::set_mesh_info( const std::vector<mesh::Entity *> &mesh_entities,
                          const VectorField * nodal_coord_ref,
                          const ScalarField * elem_weight_ref)
{
  MeshInfo mesh_info;

  /* Keep track of the total number of elements. */
  m_total_number_entities_ = mesh_entities.size();

  mesh_info.mesh_entities = mesh_entities;
  mesh_info.nodal_coord_ref = nodal_coord_ref;
  mesh_info.elem_weight_ref = elem_weight_ref;

  /** Default destination for an entity is the processor
      that already owns the entity, which is this processor.
      The length of the dest_proc_ids vector is the same
      length as the mesh_entities vector.
  */
  mesh_info.dest_proc_ids.assign(mesh_entities.size(), stk::parallel_machine_rank(comm_));

  m_mesh_information_ = mesh_info;
}

void Zoltan::init( const vector< pair<std::string,std::string> >
                   &dynamicLoadRebalancingParameters ) {
  if (0==static_zoltan_version()) {
    const double v = init_zoltan_library();
    static_zoltan_version(v);
  }

  /**
   * Create a zoltanID entity
   */

  m_zoltan_id_ = Zoltan_Create( comm_ );
  if ( m_zoltan_id_ == NULL ) {
    throw runtime_error ("(FATAL ERROR) Zoltan_Create() returned NULL");
  }

  /**
   * Set up dynamic load rebalancing
   */

  vector<pair<std::string,std::string> >::const_iterator
    P  =  dynamicLoadRebalancingParameters.begin(),
    PE =  dynamicLoadRebalancingParameters.end();

  for ( ; PE != P ; P++ ) {

    char * label = const_cast<char*>( P->first.c_str() ) ;
    char * value = const_cast<char*>( P->second.c_str() ) ;

    if (ZOLTAN_OK != (Zoltan_Set_Param(m_zoltan_id_,label,value)))
    {
      throw runtime_error(": FATAL ERROR returned from Zoltan_Set_Param ");
    }
  }

  /**
   * Register the Zoltan/SIERRA "call-back" (querry) functions.
   */
  if ( ZOLTAN_OK != register_callbacks() )
    throw runtime_error ("zoltan->Register_Callbacks error. ");

#if STK_GEOMDECOMP_DEBUG>=2
  {
    debug_print_decomp_export( zoltan, m_zoltan_id_ );
  }
#endif
  return;

}

double Zoltan::init_zoltan_library () {
  float version = 0.0;
  if ( Zoltan_Initialize( 0, NULL, &version) != ZOLTAN_OK )
    throw std::runtime_error("Return code from Zoltan_Initialize() != ZOLTAN_OK ");

  static_zoltan_version(version);
  std::ostringstream s;
  s << version;

  //sierra::ProductRegistry::instance().addTPL("Zoltan", s.str());

  return version;
}



//: ==========
//: Destructor
//: ==========

Zoltan::~Zoltan()
{
  if ( m_zoltan_id_ != NULL ) {
    Zoltan_Destroy( &m_zoltan_id_ );
  }
  m_zoltan_id_ = NULL ;
  if(Name_Conversion) {
    delete Name_Conversion;
    Name_Conversion = NULL;
  }
  if(Value_Conversion) {
    delete Value_Conversion;
    Value_Conversion = NULL;
  }
}

void
Zoltan::reset_dest_proc_data()
{
  const int  proc = 0; //Env::parallel_rank();
  const unsigned size = m_mesh_information_.mesh_entities.size();
  m_mesh_information_.dest_proc_ids.assign(size, proc);
}

void
Zoltan::set_destination_proc(const unsigned moid,
                             const unsigned proc )
{
  m_mesh_information_.dest_proc_ids[ moid ] = proc;
}

unsigned
Zoltan::destination_proc(const unsigned moid) const
{
  return m_mesh_information_.dest_proc_ids[ moid ];
}

bool
Zoltan::find_mesh_entity(const mesh::Entity * entity, unsigned & moid) const
{
  unsigned len = m_mesh_information_.mesh_entities.size();
  for(moid = 0; moid < len; ++moid)
  {
    if(m_mesh_information_.mesh_entities[moid] == entity) return true;
  }
  return false;
}

int
Zoltan::get_new_partition(stk::mesh::EntityProcVec &rebal_spec)
{
  const unsigned entity_iter_len = m_mesh_information_.mesh_entities.size();
  for (unsigned entity_iter =0; entity_iter != entity_iter_len; ++entity_iter) {
    mesh::Entity * mesh_ent = mesh_entity(entity_iter);
    int proc = destination_proc(entity_iter);
    mesh::EntityProc et(mesh_ent, proc);
    rebal_spec.push_back(et);
  }
  return 0;
}

const VectorField *
Zoltan::entity_coord_ref() const
{
   return m_mesh_information_.nodal_coord_ref;
}

mesh::Entity *
Zoltan::mesh_entity(const unsigned moid ) const
{
  return m_mesh_information_.mesh_entities[ moid ];
}

double
Zoltan::entity_weight(const unsigned moid ) const
{
  double mo_weight = 1.0;
  if (entity_weight_ref()) {
    mo_weight = * static_cast<double *>
      ( mesh::field_data (*entity_weight_ref(), *m_mesh_information_.mesh_entities[ moid ]));
  }
  return mo_weight;
}

const ScalarField *
Zoltan::entity_weight_ref() const
{
  return m_mesh_information_.elem_weight_ref;
}

unsigned
Zoltan::num_moid() const
{
  return m_mesh_information_.mesh_entities.size() ;
}

const std::string &Zoltan::parameter_entry_name() const
{
  return m_parameter_entry_name_;
}


//: Load Balance calls or Load "Re-partitioning" calls

int Zoltan::register_callbacks()
{
  /*
   * Register the Zoltan/SIERRA "call-back" (querry) functions.
   * Use ONLY THE BARE ESSENTIALS for decompositions:
   *
   *    Zoltan_Set_Num_Obj_Fn      
   *    Zoltan_Set_Obj_List_Fn     
   *    Zoltan_Set_Num_Geom_Fn
   *    Zoltan_Set_Geom_Fn
   */

  /*
   * flag for what data is to be registered with the zoltan callback
   * functions in combination with the "static_cast<CLASS>(data)"
   * statement used in the zoltan interface routine,
   * Fmwk_Zoltaninterface.C
   *
   */

  if ( Zoltan_Set_Num_Obj_Fn( m_zoltan_id_,
                              Callback_Num_Elements,
                              this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Obj_Fn using Callback_Num_Elements failed to register");
  }
  if ( Zoltan_Set_Obj_List_Fn( m_zoltan_id_, Callback_Element_List,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Obj_List_Fn using Callback_Element_List");
  }
  if ( Zoltan_Set_Num_Geom_Fn( m_zoltan_id_, Callback_Num_Dimensions,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Geom_Fn using Callback_Num_Dimensions");
  }
  if ( Zoltan_Set_Geom_Fn( m_zoltan_id_, Callback_Centroid_Coord,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Geom_Fn using Callback_Centroid_Coord");
  }
  if ( Zoltan_Set_Num_Edges_Fn( m_zoltan_id_, Callback_Num_Edges,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Num_Edges_Fn using Callback_Num_Edges");
  }
  if ( Zoltan_Set_Edge_List_Fn( m_zoltan_id_, Callback_Edge_List,this )
       != ZOLTAN_OK ) {
    throw std::runtime_error("Zoltan_Set_Edge_List_Fn using Callback_Edge_List");
  }

  return 0;

}


int  Zoltan::evaluate( int    print_stats,
                       int*   nentity,
                       double*  entity_wgt,
                       int*   ncuts,
                       double*  cut_wgt,
                       int*   nboundary,
                       int*   nadj      )
{
  int ierr        = 0;

  ZOLTAN_BALANCE_EVAL eval  = {0};
  ZOLTAN_GRAPH_EVAL   graph = {{0}};
  if (Zoltan_LB_Eval_Balance( m_zoltan_id_, print_stats, &eval)) ierr = 1;
  if (Zoltan_LB_Eval_Graph( m_zoltan_id_, print_stats, &graph) ) ierr = 1;
  *nentity      = (int)eval.nobj[0];
  *entity_wgt   =      eval.obj_wgt[0];
  *ncuts        = (int)graph.cuts[0];
  *cut_wgt      =      graph.cut_wgt[0];
  *nboundary    = (int)graph.num_boundary[0];
  *nadj         = (int)graph.nnborparts[0];

  return ierr;

}

void Zoltan::determine_new_partition (bool &RebalancingNeeded)
{
  //: Transfer export global ID lists and owning processors
  //: to SIERRA Framework's data structures

  /**
   * Perform the dynamic load rebalancing.  This just determines
   * the redistribution.  It does not actually redistribute
   * the entities.
   */

  int length_gid, length_lid;

  int           new_decomp;
  int           num_imported;
  ZOLTAN_ID_PTR import_gids;
  ZOLTAN_ID_PTR import_lids;
  int          *import_procs=NULL;
  int           num_exported;
  ZOLTAN_ID_PTR export_gids;
  ZOLTAN_ID_PTR export_lids;
  int          *export_procs=NULL;



      /** Zoltan_LB_Balance():
       *
       * lb               Pointer to the load-balancing structure, created by
       *                  Zoltan_Create, to be used in this invocation
       *                  of the load-balancing routine.
       * new_decomp       Set to 1 or .TRUE. if the decomposition was
       *                  changed by the load-balancing method; 0 or
       *                  .FALSE. otherwise.
       * num_gid_entries  Upon return, the number of array entries used to
       *                  describe a single global ID.  This value is the
       *                  maximum value over all processors of the parameter
       *                  NUM_GID_ENTRIES.
       * num_lid_entries  Upon return, the number of array entries used to
       *                  describe a single local ID.  This value is the
       *                  maximum value over all processors of the parameter
       *                  NUM_LID_ENTRIES.
       * num_imported     Upon return, the number of entities that are now
       *                  assigned to this processor that were assigned to
       *                  other processors in the old decomposition (i.e.,
       *                  the number of entities to be imported to this
       *                  processor). If the value returned is -1, no import
       *                  information has been returned and all import arrays
       *                  below are NULL (see the RETURN_LISTS parameter for
       *                  more information).
       * import_gids      Upon return, an array of num_import global IDs of
       *                  entities to be imported to this processor.
       *                  (size = num_import * num_gid_entries)
       * import_lids      Upon return, an array of num_import local IDs of
       *                  entities to be imported to this processor.
       *                  (size = num_import * num_lid_entries)
       * import_procs     Upon return, an array of size num_import listing
       *                  the processor IDs of the processors that owned the
       *                  imported entities in the previous decomposition
       *                  (i.e., the source processors).
       * num_exported     Upon return, the number of entities that were
       *                  assigned to this processor in the previous
       *                  decomposition that are now assigned to other
       *                  processors (i.e., the number of entities that must
       *                  be exported from this processor to other processors)
       *                  If the value returned is -1, no export information
       *                  has been returned and all export arrays below are
       *                  NULL (see the RETURN_LISTS parameter for more
       *                  information).
       * export_gids      Upon return, an array of num_export global IDs of
       *                  entities to be exported from this processor.
       *                  (size = num_export * num_gid_entries)
       * export_lids      Upon return, an array of num_export local IDs of
       *                  entities to be exported from this processor.
       *                  (size = num_export * num_lid_entries)
       * export_procs     Upon return, an array of size num_export listing
       *                  the processor IDs of processors that will own the
       *                  exported entities in the new decomposition (i.e.,
       *                  the destination processors).
       */

  int status = Zoltan_LB_Balance( m_zoltan_id_,        &new_decomp,
                                  &length_gid     , &length_lid,
                                  &num_imported,    &import_gids,
                                  &import_lids,     &import_procs,
                                  &num_exported,    &export_gids,
                                  &export_lids,     &export_procs );
  if (status != ZOLTAN_OK) {
    throw std::runtime_error("Zoltan_Balance() returned error code " + status);
  }

  //: Initialize destination processor IDs (dest_proc_ids)
  reset_dest_proc_data();

  int actual_exported = 0;
  if ( new_decomp && ( num_exported != -1 ) ) {
    const unsigned parallel_rank = parallel_machine_rank(comm_);
    /* New Decomposition was generated */
    for (int j=0; j < num_exported; ++j ) {

      //: Get exported region, local, global, and processor ids
      //const unsigned rid = export_gids[ j*::num_gid_entries() ];  // Region ID variable
      //if (!rid) throw runtime_error ("Region ID variable should be non-zero.");
      const unsigned lid = export_lids[ j*::num_lid_entries() ];  // Local  ID variable
      const unsigned pid = export_procs[ j ];                     // Exported Processor ID (i.e., MPI "rank" )

      if (parallel_rank != pid) {
        ++actual_exported;
        set_destination_proc(lid, pid);
      }
    }
  }

  RebalancingNeeded = 0 ;
  if (new_decomp) {
    int rebalneeded=0;
    stk::all_reduce_sum(comm_, &actual_exported, &rebalneeded, 1);
    if (rebalneeded)  RebalancingNeeded = 1;
  }

  /**
   * Clean up after zoltan
   */
  if ( ZOLTAN_OK !=
       Zoltan_LB_Free_Data( &import_gids, &import_lids, &import_procs,
                            &export_gids, &export_lids, &export_procs )) {
    throw runtime_error (" FATAL ERROR in Zoltan_LB_Free_Data.");
  }

}

void Zoltan::convert_names_and_values(const Parameters &from, Parameters &to)
{
  /* First time through, fill the conversion tables. */
  if (!Name_Conversion) {
    Name_Conversion = new Parameters;
    fill_name_conversion (*Name_Conversion);
  }
  if (!Value_Conversion) {
    Value_Conversion = new Parameters;
    fill_value_conversion(*Value_Conversion);
  }

  // NOTE: "LOAD BALANCING METHOD" default
  // is also hard coded in fill_default_values();
  std::string algorithm;
  const std::string keyname("LOAD BALANCING METHOD");
  if( from.isParameter(keyname) )
    algorithm = from.get<std::string>(keyname);
  else
    algorithm = "0";

  const Parameters & General = Name_Conversion->sublist("General");
  const Parameters & Algorithm = Name_Conversion->sublist(algorithm);

  Parameters::ConstIterator
    from_iter  = from.begin(),
    from_end   = from.end();

  /* Iterate over all of the input parameters to find proper
     Zoltan names and Zoltan parameter values. */
  for (; from_iter != from_end; ++from_iter) {
    std::string
      from_name = from.name(from_iter),
      to_name;

    /* Check to see if this is a general parameter name
       and if not, check if it is an algorithm specific name. */

    if (General.isParameter(from_name)) {
      to_name = General.get<std::string>(from_name);
    } else {
      to_name = Algorithm.get<std::string>(from_name);
    }

    /* Now convert the parameter value to the correct form.
       Only a couple of parameters have parameter conversion.
       The ones converted are nested in Value_Conversion.
    */
    std::string to_value = Teuchos::getValue<string>(from.entry(from_iter));
    //if (Value_Conversion->isParameter(from_name)) to_value = Value_Conversion->get<std::string>(to_value);
    if (Value_Conversion->isParameter(from_name)) to_value = Value_Conversion->sublist(from_name).get<std::string>(to_value);
    if (!to_name.empty()) to.set(to_name, to_value);
  }
}

void Zoltan::merge_default_values(const Parameters &from,
                                        Parameters &to)
{
  Parameters default_values;
  fill_default_values(default_values);
  to.setParameters(default_values);
  to.setParameters(from);
}
