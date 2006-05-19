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
    NumVertices_(0),
  NumDimensions_(0),
  NumComponents_(0)
  {}

  //! Copy constructor.
  Element(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumDimensions(rhs.getNumDimensions());
    setNumComponents(rhs.getNumComponents());
    for (int i = 0; i < rhs.getNumComponents(); ++i)
      setComponent(i, rhs.getComponent(i));
    setLabel(rhs.getLabel());
  }

  Element& operator=(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumDimensions(rhs.getNumDimensions());
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
    return(NumVertices_);
  }

  int getNumDimensions() const
  {
    return(NumDimensions_);
  }

  int getNumComponents() const 
  {
    return(NumComponents_);
  }

  const Element& getComponent(const int which) const
  {
    return(Components_[which]);
  }

  // @}
  // @{ \name Set methods

  void setNumVertices(const int NumVertices)
  {
    NumVertices_ = NumVertices;
  }

  void setNumDimensions(const int NumDimensions)
  {
    NumDimensions_ = NumDimensions;
  }

  void setNumComponents(const int NumComponents)
  {
    NumComponents_ = NumComponents;
    Components_.resize(NumComponents);
  }

  void setComponent(const int which, const Element& what)
  {
    Components_[which] = what;
  }

  virtual void Print(ostream & os) const
  {
    cout << endl;
    cout << "Number of vertices = " << getNumVertices() << endl;
    cout << "Number of dimensions = " << getNumDimensions() << endl;
    cout << "Number of components = " << getNumComponents() << endl;
    for (int i = 0; i < getNumComponents(); ++i)
      cout << "Label of component " << i << " = " <<
        getComponent(i).getLabel() << endl;
    cout << endl;
  }

  // @}
private:
  // @{ \name Private data and methods

  int NumVertices_;
  int NumDimensions_;
  int NumComponents_;
  vector<Element> Components_;

  // @}
}; // class Element

} // namespace grid
} // namespace phx

#endif // PHX_BASE_ELEMENT_H
