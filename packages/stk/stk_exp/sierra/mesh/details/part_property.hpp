#ifndef SIERRA_SIERRA_MESH_MESH_PART_HPP
#define SIERRA_SIERRA_MESH_MESH_PART_HPP


#include <sierra/mesh/details/part_key.hpp>

#include <boost/any.hpp>

#include <string>
#include <typeinfo>
#include <stdexcept>
#include <vector>
#include <cstdlib>

namespace sierra {
namespace mesh {
namespace details {

class part_property
{
  public:

    part_property()
      : m_name("INVALID_PART"), m_key(), m_properties()
    {}

    part_property(const std::string & arg_name, part_key arg_key )
      : m_name(arg_name), m_key(arg_key), m_properties()
    {}

    template<class T>
    part_property( const std::string & arg_name, part_key arg_key, const T& arg_prop)
      : m_name(arg_name), m_key(arg_key), m_properties(1, boost::any(arg_prop))
    {
    }

    ~part_property()
    {
    }

    bool operator == ( const part_property & rhs) const
    {
      return m_key == rhs.m_key;
    }

    bool operator != ( const part_property & rhs) const
    {
      return !(*this == rhs);
    }

    const std::string & name() const { return m_name; }

    part_key key() const { return m_key; }


    template<class T>
    void add_property(const T& prop)
    {
      for(size_t i=0; i<m_properties.size(); ++i) {
        const T* value = boost::any_cast<T>(&m_properties[i]);
        if (value != NULL) {
          m_properties[i] = prop;
          return;
        }
      }
      m_properties.push_back(boost::any(prop));
    }

    template<class T>
    bool has_property() const
    {
      for(size_t i=0; i<m_properties.size(); ++i) {
        if (boost::any_cast<T>(&m_properties[i]) != NULL) return true;
      }
      return false;
    }

    template<class T>
    const T& get_property() const
    {
      const T* value = NULL;
      for(size_t i=0; i<m_properties.size(); ++i) {
        value = boost::any_cast<T>(&m_properties[i]);
        if (value != NULL) break;
      }

      if (value == NULL) throw std::runtime_error("sierra::mesh::details::part_property::get_property failed to find property");
      return *value;
    }

  private:
    std::string        m_name;
    part_key           m_key;
    std::vector<boost::any>   m_properties;

};

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_MESH_PART_HPP
