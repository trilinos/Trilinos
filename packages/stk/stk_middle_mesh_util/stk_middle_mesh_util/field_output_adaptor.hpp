#ifndef STK_MIDDLE_MESH_UTIL_FIELD_OUTPUT_ADAPTOR_H
#define STK_MIDDLE_MESH_UTIL_FIELD_OUTPUT_ADAPTOR_H

#include "stk_middle_mesh/field.hpp"
#include <string>
#include <memory>


#include <stk_io/IossBridge.hpp>

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

class FieldOutputAdaptor
{
  public: 
    virtual ~FieldOutputAdaptor() = default;

    virtual std::string get_name() const = 0;

    virtual stk::middle_mesh::mesh::FieldShape get_field_shape() const = 0;

    virtual int get_num_comp() const = 0;

    virtual void get_data(mesh::MeshEntityPtr entity, int node, std::vector<double>& vals) const = 0;

    virtual stk::io::FieldOutputType get_field_output_type() const { return stk::io::FieldOutputType::SCALAR; }
};

using FieldOutputAdaptorPtr = std::shared_ptr<FieldOutputAdaptor>;

class FieldOutputAdaptorDouble : public FieldOutputAdaptor
{
  public:
    FieldOutputAdaptorDouble(stk::middle_mesh::mesh::FieldPtr<double> field, const std::string& name) :
      m_field(field),
      m_name(name)
    {}

    std::string get_name() const override { return m_name; }

    stk::middle_mesh::mesh::FieldShape get_field_shape() const override {return m_field->get_field_shape(); }

    int get_num_comp() const override { return m_field->get_num_comp(); }

    void get_data(mesh::MeshEntityPtr entity, int node, std::vector<double>& vals) const override
    {
      vals.resize(m_field->get_num_comp());

      auto& field = *m_field;
      for (int i=0; i < field.get_num_comp(); ++i)
        vals[i] = field(entity, node, i);
    }

  private:
    stk::middle_mesh::mesh::FieldPtr<double> m_field;
    std::string m_name;
};

}
}
}
}

#endif