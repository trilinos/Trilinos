#include <iomanip>
#include <sstream>
#include <iostream>

#include <stk_sddm/DataTypeTraits.hpp>

namespace stk {
namespace sddm {

const char *
DataTypeTraits<int>::name() 
{  
  return "Integer";
}

std::ostream &
DataTypeTraits<int>::dump(
  std::ostream &        os,
  const int &           t)
{
  os << t;
  return os;
}

std::ostream &
DataTypeTraits<int>::xml(
  std::ostream &        os,
  const int &           t)
{
  os << "<Integer>" << t << "</Integer>";
  return os;
}

std::istream &
DataTypeTraits<int>::load(
  std::istream &        is,
  int &                 t)
{
  is >> t;
  return is;
}

const char *
DataTypeTraits<double> ::name() 
{
  return "Double";
}

std::ostream &
DataTypeTraits<double>::dump(
  std::ostream &        os,
  const double &        t)
{
  std::streamsize p = os.precision();
  os.precision(16);
  os << t;
  os.precision(p);
  
  return os;
}

std::istream &
DataTypeTraits<double>::load(
  std::istream &        is,
  double &              t)
{
  is >> t;
  return is;
}

namespace {


std::ostream &
encode(
  std::ostream &        os,
  const std::string &   s)
{
  for (std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
    switch (*it) {    
    case '\n': os << "\n"; break;
    case '\t': os << "\t"; break;
    case '\v': os << "\v"; break;
    case '\b': os << "\b"; break;
    case '\r': os << "\r"; break;
    case '\f': os << "\f"; break;
    case '\a': os << "\a"; break;
    case '\?': os << "\?"; break;
    case '\\': os << "\\"; break;
    case '\'': os << "\'"; break;
    case '\"': os << "\""; break;
    default:
      if (*it < ' ' || *it == '\127') {
        char fill_char = os.fill();
        std::ios_base::fmtflags fmt_flags = os.flags();
        
        os << "\\";
        os.width(3);
        os.fill('0');
        os.setf(std::ios_base::oct, std::ios_base::basefield);
        os << (unsigned int) (*it);
        os.fill(fill_char);
        os.flags(fmt_flags);
      }
      else
        os << *it;
      break;
    }
  }
  return os;
}

std::istream &
decode(
  std::istream &        is,
  std::string &         s) 
{

  char c = is.get();
  
  if (c != '\"') {
    is.putback(c);
    is >> s;
  }
  else {
    std::ostringstream oss;
    
    for (;;) {
      c = is.get();
      if (c == '\"')
        break;
      else if (c == '\\') {
        char d = is.get();
        switch (d) {
        case 'n': oss << '\n'; break;
        case 't': oss << '\t'; break;
        case 'v': oss << '\v'; break;
        case 'b': oss << '\b'; break;
        case 'r': oss << '\r'; break;
        case 'f': oss << '\f'; break;
        case 'a': oss << '\a'; break;
        case '?': oss << '\?'; break;
        case '\\': oss << '\\'; break;
        case '\'': oss << '\''; break;
        case '\"': oss << '\"'; break;
        case '0': 
          {
            int i = (is.get() - '0')*0100 + (is.get() - '0')*0x10 + (is.get() - '0');
            oss << (char) i;
          }
        }
      }
    }
    is.get();

    s = oss.str();
  }
  return is;
}

} // namespace <unnamed>


const char *
DataTypeTraits<std::string>::name() 
{
  return "String";
}

std::ostream &
DataTypeTraits<std::string>::dump(
  std::ostream &        os,
  const std::string &   t)
{
  encode(os, t);
  
  return os;
}


std::istream &
DataTypeTraits<std::string>::load(
  std::istream &        is,
  std::string &         t)
{
  return decode(is, t);
}


const char *
DataTypeTraits<std::vector<int> >::name() 
{
  return "ListOfInteger";
}

// std::ostream &
// DataTypeTraits<std::vector<int> >::dump(
//   std::ostream &                os,
//   const std::vector<int> &      t)
// {
//   for (std::vector<int>::const_iterator it = t.begin(); it != t.end(); ++it) {
//     if (it != t.begin())
//       os << " ";
//     DataTypeTraits<int>::dump(os, (*it));
//   }
  
//   return os;
// }

// std::istream &
// DataTypeTraits<std::vector<int> >::load(
//   std::istream &        is,
//   std::vector<int> &    t)
// {
//   t.clear();
  
//   while (!is.eof()) {
//     int v;
//     DataTypeTraits<int>::load(is, v);
//     t.push_back(v);
//   }
  
//   return is;
// }


const char *
DataTypeTraits<std::vector<double> >::name() 
{
  return "ListOfDouble";
}

// std::ostream &
// DataTypeTraits<std::vector<double> >::dump(
//   std::ostream &                os,
//   const std::vector<double>  &  t)
// {
//   for (std::vector<double>::const_iterator it = t.begin(); it != t.end(); ++it) {
//     if (it != t.begin())
//       os << " ";
//     DataTypeTraits<double>::dump(os, (*it));
//   }
  
//   return os;
// }

// std::istream &
// DataTypeTraits<std::vector<double> >::load(
//   std::istream &        is,
//   std::vector<double> & t)
// {
//   t.clear();
  
//   while (!is.eof()) {
//     double v;
//     DataTypeTraits<double>::load(is, v);
//     t.push_back(v);
//   }
  
//   return is;
// }


const char *
DataTypeTraits<std::vector<std::string> >::name() 
{
  return "ListOfString";
}


// std::ostream &
// DataTypeTraits<std::vector<std::string> >::dump(
//   std::ostream &                        os,
//   const std::vector<std::string>  &     t)
// {
//   for (std::vector<std::string>::const_iterator it = t.begin(); it != t.end(); ++it) {
//     if (it != t.begin())
//       os << " ";
//     os << "\"";
//     DataTypeTraits<std::string>::dump(os, (*it));
//     os << "\"";
//   }
  
//   return os;
// }

// std::istream &
// DataTypeTraits<std::vector<std::string> >::load(
//   std::istream &                is,
//   std::vector<std::string> &    t)
// {
//   t.clear();
  
//   while (!is.eof()) {
//     std::string v;
//     DataTypeTraits<std::string>::load(is, v);
//     t.push_back(v);
//   }
  
//   return is;
// }

} // namespace sddm
} // namespace stk
