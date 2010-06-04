/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_ElementTopology.h>
#include <Ioss_Utils.h>
#include <Ioss_Super.h>

#include <cstring>
#include <assert.h>


// ========================================================================
Ioss::ElementTopology::ElementTopology(const std::string &type, const std::string &master_elem_name)
  : name_(type), masterElementName_(master_elem_name)
{
  registry()->insert(Ioss::ETM_VP(type, this));
  alias(name_, masterElementName_);
}

void Ioss::ElementTopology::alias(const std::string& base, const std::string& syn)
{
  registry()->insert(Ioss::ETM_VP(syn, factory(base)));
}

Ioss::ElementTopologyMap* Ioss::ElementTopology::registry()
{
  static Ioss::ElementTopologyMap registry_ ;
  return &registry_;
}

Ioss::ElementTopology::~ElementTopology() {}

bool Ioss::ElementTopology::edges_similar() const {return true;}
bool Ioss::ElementTopology::faces_similar() const {return true;}

Ioss::ElementTopology* Ioss::ElementTopology::factory(const std::string& type, bool ok_to_fail)
{
  Ioss::ElementTopology* inst = NULL;
  Ioss::ElementTopologyMap::iterator iter = registry()->find(type);
  if (iter == registry()->end() && std::strncmp(type.c_str(), "super", 5) == 0) {
    // A super element can have a varying number of nodes.  Create
    // an IO element type for this super element. The node count
    // should be encoded in the 'type' as 'super42' for a 42-node
    // superelement.
    
    Ioss::Super::make_super(type);
    iter = registry()->find(type);
  }

  if (iter == registry()->end()) {
    if (!ok_to_fail)
      IOSS_WARNING << "WARNING: The topology type '" << type << "' is not supported.\n";
  } else {
    inst = (*iter).second;
  }
  return inst;
}

Ioss::ElementTopology* Ioss::ElementTopology::factory(unsigned int unique_id)
{
  // Given a unique id obtained from 'get_unique_id', return the
  // topology type that it refers to...
  Ioss::ElementTopologyMap::const_iterator I;
  for (I = registry()->begin(); I != registry()->end(); ++I) {
    if (Ioss::Utils::hash((*I).second->name()) == unique_id)
      return (*I).second;
  }
  return NULL;
}

unsigned int Ioss::ElementTopology::get_unique_id(const std::string& type)
{
  // Return a unique integer id corresponding to this topology type.
  // Basically used to simplify some parallel calculations so they can
  // deal with int instead of strings...
  if (type == "unknown")
    return 0;
  unsigned int hash_val = 0;
  Ioss::ElementTopologyMap::iterator iter = registry()->find(type);
  if (iter == registry()->end()) {
    IOSS_WARNING << "WARNING: The topology type '" << type
		 << "' is not supported.\n";
  } else {
    Ioss::ElementTopology* inst = (*iter).second;
    hash_val = Ioss::Utils::hash(inst->name());
  }
  return hash_val;
}

int Ioss::ElementTopology::describe(NameList *names)
{
  int count = 0;
  Ioss::ElementTopologyMap::const_iterator I;
  for (I = registry()->begin(); I != registry()->end(); ++I) {
    names->push_back((*I).first);
    count++;
  }
  return count;
}

Ioss::IntVector Ioss::ElementTopology::face_edge_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= number_faces());

  int nface_edge = number_edges_face(face_number);
  Ioss::IntVector fcon(nface_edge);

  // This works for 2D elements, 3D elements override
  for (int i=0; i < nface_edge; i++)
    fcon[i] = i;

  return fcon;
}

Ioss::IntVector Ioss::ElementTopology::element_edge_connectivity() const
{
  int nedge = number_edges();
  Ioss::IntVector econ(nedge);
  for (int i=0; i < nedge; i++)
    econ[i] = i;

  return econ;
}

bool Ioss::ElementTopology::is_alias(const std::string &my_alias) const
{
  Ioss::ElementTopologyMap::iterator iter = registry()->find(my_alias);
  if (iter == registry()->end()) {
    return false;
  } else {
    return this == (*iter).second;
  }
}

bool Ioss::ElementTopology::is_element() const
{
  // NOTE: This is overridden in some derived classes.
  // The definition here is the default if not overridden.
  if (spatial_dimension() == parametric_dimension())
    return true;
  else
    return false;
}

int Ioss::ElementTopology::number_boundaries() const
{
  if (parametric_dimension() == 3 && spatial_dimension() == 3)
    return number_faces();

  if (parametric_dimension() == 2 && spatial_dimension() == 2)
    return number_edges();

  if (is_element()) {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return number_faces();
    } else if (parametric_dimension() == 1) {
      return number_edges();
    }
  } else {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return number_edges();
    }
  }    
  return 0;
}

Ioss::IntVector Ioss::ElementTopology::boundary_connectivity(int edge_number) const
{
  if (parametric_dimension() == 3 && spatial_dimension() == 3)
    return face_connectivity(edge_number);

  if (parametric_dimension() == 2 && spatial_dimension() == 2)
    return edge_connectivity(edge_number);

  if (is_element()) {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return face_connectivity(edge_number);
    } else if (parametric_dimension() == 1) {
      return edge_connectivity(edge_number);
    }
  } else {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return edge_connectivity(edge_number);
    }
  }    
  return Ioss::IntVector();
}

Ioss::ElementTopology* Ioss::ElementTopology::boundary_type(int face_number) const
{
  if (parametric_dimension() == 3 && spatial_dimension() == 3)
    return face_type(face_number);

  if (parametric_dimension() == 2 && spatial_dimension() == 2)
    return edge_type(face_number);

  if (is_element()) {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return face_type(face_number);
    } else if (parametric_dimension() == 1) {
      return edge_type(face_number);
    }
  } else {
    if (parametric_dimension() == 2) {
      assert(spatial_dimension() == 3);
      return edge_type(face_number);
    }
  }    
  return NULL;
}
