#ifndef PHX_GRID_ELEMENT_H
#define PHX_GRID_ELEMENT_H

#include "Galeri_ConfigDefs.h"
#include "phx_core_Object.h"
#include <vector>

namespace phx {
namespace grid {

class Element : public phx::core::Object
{
public: 
  // @{ \name constructor and destructor
  //! Default constructor.
  Element() :
    numVertices_(0),
    numComponents_(0)
  {}

  //! Copy constructor.
  Element(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumComponents(rhs.getNumComponents());
    for (int i = 0; i < rhs.getNumComponents(); ++i)
      setComponent(i, rhs.getComponent(i));
    setLabel(rhs.getLabel());
  }

  Element& operator=(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumComponents(rhs.getNumComponents());
    for (int i = 0; i < rhs.getNumComponents(); ++i)
      setComponent(i, rhs.getComponent(i));
    setLabel(rhs.getLabel());

    return(*this);
  }

  //! destructor.
  ~Element() {}

  // @}
  // @{ \name Get methods

  int getNumVertices() const
  {
    return(numVertices_);
  }

  int getNumComponents() const 
  {
    return(numComponents_);
  }

  const Element& getComponent(const int which) const
  {
    return(components_[which]);
  }

  // @}
  // @{ \name Set methods

  void setNumVertices(const int numVertices)
  {
    numVertices_ = numVertices;
  }

  void setNumComponents(const int numComponents)
  {
    numComponents_ = numComponents;
    components_.resize(numComponents);
  }

  void setComponent(const int which, const Element& what)
  {
    components_[which] = what;
  }

  virtual void print(ostream & os) const
  {
    cout << endl;
    cout << "Number of vertices = " << getNumVertices() << endl;
    cout << "Number of components = " << getNumComponents() << endl;
    for (int i = 0; i < getNumComponents(); ++i)
      cout << "Label of component " << i << " = " <<
        getComponent(i).getLabel() << endl;
    cout << endl;
  }

  // @}
private:
  // @{ \name Private data and methods

  int numVertices_;
  int numComponents_;
  vector<Element> components_;

  // @}
}; // class Element

} // namespace grid
} // namespace phx
#endif // PHX_BASE_ELEMENT_H
