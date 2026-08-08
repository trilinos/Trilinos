#ifndef PANZER_STK_MESH_ACCESSOR_HPP
#define PANZER_STK_MESH_ACCESSOR_HPP

#include "Teuchos_Assert.hpp"
#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

  /** Base class for closure model factories to provide access to the STK_Interface. 

      Applications can derive their closure model factories from this
      class. The model evalautor factory will create and set the mesh
      on any factories that derive from this object before calling the
      buildObjects().
   */
  class STKMeshAccessor {
    Teuchos::RCP<panzer_stk::STK_Interface> mesh_;
  public:
    /// \brief Returns the mesh set by setMesh(). Asserts that a mesh has been set.
    Teuchos::RCP<panzer_stk::STK_Interface> getMesh() const
    {
      TEUCHOS_ASSERT(Teuchos::nonnull(mesh_));
      return mesh_;
    }

    /// \brief Sets the mesh to be returned by getMesh().
    void setMesh(const Teuchos::RCP<STK_Interface>& mesh){mesh_ = mesh;}
  };

}

#endif
