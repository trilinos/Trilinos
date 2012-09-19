#ifndef SDDM_DATATYPES_HPP
#define SDDM_DATATYPES_HPP

#include <iosfwd>
#include <vector>

namespace stk {
namespace sddm {

template<class T>
struct DataTypeTraits;


template<>
struct DataTypeTraits<int> 
{
  typedef int Type;

  static const char *name();
  static std::ostream &dump(std::ostream &os, const Type &t);
  static std::ostream &xml(std::ostream &os, const Type &t);
  static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<float> 
{
  typedef float Type;

  static const char *name();
  static std::ostream &dump(std::ostream &os, const Type &t);
  static std::ostream &xml(std::ostream &os, const Type &t);
  static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<double> 
{
  typedef double Type;

  static const char *name();
  static std::ostream &dump(std::ostream &os, const Type &t);
  static std::ostream &xml(std::ostream &os, const Type &t);
  static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<std::string> 
{
  typedef std::string Type;

  static const char *name();
  static std::ostream &dump(std::ostream &os, const Type &t);
  static std::ostream &xml(std::ostream &os, const Type &t);
  static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<std::vector<int> > 
{
  typedef std::vector<int> Type;

  static const char *name();
  // static std::ostream &dump(std::ostream &os, const Type &t);
  // static std::ostream &xml(std::ostream &os, const Type &t);
  // static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<std::vector<double> > 
{
  typedef std::vector<double> Type;

  static const char *name();
  // static std::ostream &dump(std::ostream &os, const Type &t);
  // static std::ostream &xml(std::ostream &os, const Type &t);
  // static std::istream &load(std::istream &is, Type &t);
};

template<>
struct DataTypeTraits<std::vector<std::string> > 
{
  typedef std::vector<std::string> Type;

  static const char *name();
  // static std::ostream &dump(std::ostream &os, const Type &t);
  // static std::ostream &xml(std::ostream &os, const Type &t);
  // static std::istream &load(std::istream &is, Type &t);
};

} // namespace sddm
} // namespace stk

#endif // SDDM_DATATYPES_HPP
