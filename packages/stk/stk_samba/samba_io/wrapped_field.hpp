#ifndef SAMBA_SAMBA_IO_FIELD_WRAPPER_HPP
#define SAMBA_SAMBA_IO_FIELD_WRAPPER_HPP

#include <samba_io/ioss_field_helper.hpp>

namespace samba {
namespace io {

class wrapped_field_base
{
public:

    wrapped_field_base(const std::string &io_name) : m_io_name(io_name) { }
    virtual ~wrapped_field_base() { }

    virtual Ioss::Field::BasicType basicType() const = 0;
    virtual std::string storage(mesh mesh_arg) const = 0;
    virtual size_t valueCount() const = 0;

    std::string getName() const { return m_io_name; }

private:

    std::string m_io_name;

};


template <typename Field_T>
class wrapped_field : public wrapped_field_base
{
public:

    wrapped_field(Field_T field_arg, mesh mesh_arg, const std::string &io_name)
        : wrapped_field_base(io_name)
        , m_field(field_arg)
        , m_mesh(mesh_arg)
    {}

    virtual Ioss::Field::BasicType basicType() const { return computeBasicType(); }

    virtual std::string storage(mesh mesh_arg) const { return computeStorage(mesh_arg); }

    virtual size_t valueCount() const { return m_mesh.num_entities(m_field.restriction());}

private:

    Field_T m_field;
    mesh m_mesh;

    detail::IossFieldHelper<typename Field_T::data_type, typename Field_T::dimension_functor>
        m_configure_field;

    Ioss::Field::BasicType computeBasicType() const { return m_configure_field.basicType();}
    std::string computeStorage(mesh mesh_arg) const { return m_configure_field.storage(); }

};

} // namespace io
} // namespace samba


#endif // SAMBA_SAMBA_IO_FIELD_WRAPPER_HPP
