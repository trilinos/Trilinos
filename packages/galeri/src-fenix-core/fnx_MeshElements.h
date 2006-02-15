#ifndef FNX_ELEMENTS_H
#define FNX_ELEMENTS_H

#include "Epetra_Object.h"

namespace fnx
{
  class MeshElement : public Epetra_Object
  {
    public:
      virtual ~MeshElement() {}

      virtual int NumDimensions() const = 0;

      virtual int NumEntities(const int dim) const = 0;
  };

  class Triangle : public MeshElement
  {
    public:
      virtual Triangle() : Epetra_Object("fnx::Element::Triangle") {}

      virtual int NumDimensions() const
      {
        return(2);
      }

      virtual int NumEntities(const int dim) const
      {
        if (dim == 0)      return(3);
        else if (dim == 1) return(3);
        else               FNX_THROW("Dimension not valid");
      }
  };

  class Quad : public MeshElement
  {
    public:
      virtual Quad() : Epetra_Object("fnx::Element::Quad") {}

      virtual int NumDimensions() const
      {
        return(2);
      }

      virtual int NumEntities(const int dim) const
      {
        if (dim == 0)      return(4);
        else if (dim == 1) return(4);
        else               FNX_THROW("Dimension not valid");
      }
  };

  class Tet : public MeshElement
  {
    public:
      virtual Tet() : Epetra_Object("fnx::Element::Tet") {}

      virtual int NumDimensions() const
      {
        return(3);
      }

      virtual int NumEntities(const int dim) const
      {
        if (dim == 0)      return(4);
        else if (dim == 1) return(6);
        else if (dim == 2) return(4);
        else               FNX_THROW("Dimension not valid");
      }
  };

  class Hex : public MeshElement
  {
    public:
      virtual Hex() : Epetra_Object("fnx::Element::Hex") {}

      virtual int NumDimensions() const
      {
        return(3);
      }

      virtual int NumEntities(const int dim) const
      {
        if (dim == 0)      return(8);
        else if (dim == 1) return(12);
        else if (dim == 2) return(6);
        else               FNX_THROW("Dimension not valid");
      }
  };
  
}; // namespace fnx

#endif
