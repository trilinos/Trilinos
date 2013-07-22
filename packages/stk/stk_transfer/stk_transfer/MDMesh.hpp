
#include <boost/shared_ptr.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>


namespace stk {
namespace transfer {

template <unsigned DIM> class MDMesh {
public :
  typedef Intrepid::FieldContainer<double> MDArray;
  typedef unsigned                                           Entity;
  typedef std::vector<Entity>                                EntityVec;
  typedef unsigned                                           EntityKey;
  typedef std::set   <EntityKey>                             EntityKeySet;
  typedef search::ident::IdentProc<EntityKey, unsigned>      EntityProc;
  typedef std::vector<EntityProc>                            EntityProcVec;

  typedef search::box::SphereBoundingBox<EntityProc,float,DIM> BoundingBox;

  enum {Dimension = DIM};

  MDMesh(MDArray                     &val,
         const MDArray               &coord,
         const double                 initial_radius,
         const stk::ParallelMachine   comm);
  ~MDMesh();

  // Needed for STK Transfer
  ParallelMachine comm() const {return m_comm;}

  void bounding_boxes (std::vector<BoundingBox> &v) const;

  void update_values();

  // Needed for LinearInterpoate
  const double *coord(const EntityKey k) const;
  const double *value(const EntityKey k, const unsigned i=0) const;
        double *value(const EntityKey k, const unsigned i=0);
  unsigned      value_size(const EntityKey e, const unsigned i=0) const;
  unsigned      num_values() const;

  struct Record { virtual ~Record(){} };
  template <class T> T* database(const EntityKey k) {
    typename RecordMap::const_iterator i = m_record_map.find(k);
    if (i == m_record_map.end()) {
      RecordPtr record(new T());
      typename RecordMap::value_type v(k,record);
      i = m_record_map.insert(v).first;
    }
    T *record = dynamic_cast<T*>(i->second.get());
    ThrowRequireMsg (record,__FILE__<<":"<<__LINE__<<" Dynamic Cast failed in MDMesh::database ");
    return record;
  }
  

private :
  MDMesh (); 
  MDMesh(const MDMesh &M);
  MDMesh &operator=(const MDMesh&);

  const unsigned                         m_num_nodes;
  const unsigned                       m_spatial_dim;
  const unsigned                        m_num_values;
  const double                          m_sphere_rad;
  const MDArray                 &m_coordinates_field;
        MDArray                      &m_values_field;
  const ParallelMachine                       m_comm;

  typedef boost::shared_ptr<Record>         RecordPtr;
  typedef std::map<EntityKey,RecordPtr>     RecordMap;
  RecordMap                                 m_record_map;
};

template<unsigned DIM> MDMesh<DIM>::MDMesh(
                MDArray               &val,
          const MDArray               &coord,
          const double                 initial_radius,
          const stk::ParallelMachine   comm) :
    m_num_nodes         (coord.dimension(0)),
    m_spatial_dim       (coord.dimension(1)),
    m_num_values        (val  .dimension(1)), 
    m_sphere_rad        (initial_radius    ), 
    m_coordinates_field (coord), 
    m_values_field      (val)  ,
    m_comm              (comm)
{}

template<unsigned DIM> MDMesh<DIM>::~MDMesh(){}

template<unsigned DIM> void MDMesh<DIM>::bounding_boxes (std::vector<BoundingBox> &v) const {

  typedef typename BoundingBox::Data Data;
  typedef typename BoundingBox::Key  Key;
  const Data r=m_sphere_rad;

  v.clear();

  for (unsigned Id=0; Id!=m_num_nodes; ++Id) {

    Data center[Dimension];
    const double *c = &m_coordinates_field(Id,0);
    for (unsigned j=0; j<Dimension; ++j) center[j] = c[j];
    const Key key(Id, parallel_machine_rank(comm()));
    BoundingBox B(center, r, key);
    v.push_back(B);
  }
}

template<unsigned DIM> void MDMesh<DIM>::update_values () {}
  
template<unsigned DIM> const double *MDMesh<DIM>::coord(const EntityKey k) const {
  const double *c = &m_coordinates_field(k,0);
  return  c;
}

template<unsigned DIM> unsigned  MDMesh<DIM>::num_values() const {
 return m_num_values;
}

template<unsigned DIM> unsigned  MDMesh<DIM>::value_size(const EntityKey k, const unsigned i) const {
  return  1;
}

template<unsigned DIM> const double *MDMesh<DIM>::value(const EntityKey k, const unsigned i) const {
  const double *value = &m_values_field(k, i);
  return  value;
}

template<unsigned DIM> double *MDMesh<DIM>::value(const EntityKey k, const unsigned i) {
  double *value = &m_values_field(k, i);
  return  value;
}

}
}
