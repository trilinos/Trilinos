/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/impl_NamedValues.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk
{
namespace coupling
{
namespace impl
{

class Pack
{
public:
  Pack(stk::CommBuffer& buffer)
    : m_buffer(buffer)
  { }

  std::string name() { return "pack"; }

  template <typename T>
  void execute(Value& value)
  {
    m_buffer.pack(std::any_cast<T>(value.value));
  }

private:
  stk::CommBuffer& m_buffer;
};

void NamedValues::pack(stk::CommBuffer& buf) const
{
  buf.pack<unsigned>(m_values.size());
  for(auto&& iter : m_values) {
    const std::string& name = iter.first;
    buf.pack(name);

    Value value = iter.second;
    buf.pack<unsigned char>(static_cast<unsigned char>(value.type));

    Pack packer(buf);
    execute(packer, value);
  }
}

class Unpack
{
public:
  Unpack(stk::CommBuffer& buffer)
    : m_buffer(buffer)
  { }

  std::string name() { return "unpack"; }

  template <typename T>
  void execute(Value& value)
  {
    T val;
    m_buffer.unpack(val);
    value = Value(val);
  }

private:
  stk::CommBuffer& m_buffer;
};

void NamedValues::unpack(stk::CommBuffer& buf)
{
  unsigned numItems = 0;
  buf.unpack<unsigned>(numItems);
  for (unsigned i=0; i<numItems; ++i) {
    std::string name;
    buf.unpack(name);

    Value value;

    unsigned char c;
    buf.unpack<unsigned char>(c);
    value.type = static_cast<ValueTypes>(c);

    Unpack unpacker(buf);
    execute(unpacker, value);

    m_values.insert(std::make_pair(name, value));
  }
}

} // namespace impl
} // namespace coupling
} // namespace stk
