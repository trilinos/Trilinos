#ifndef ADR_PROBLEM_HPP
#define ADR_PROBLEM_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_ConfigDefs.h"

namespace ADR {

namespace Xpetra {

enum {
  DIR_LEFT   = 0x01,
  DIR_RIGHT  = 0x02,
  DIR_BOTTOM = 0x04,
  DIR_TOP    = 0x08,
  DIR_FRONT  = 0x10,
  DIR_BACK   = 0x20,
  DIR_ALL    = DIR_LEFT | DIR_RIGHT | DIR_BOTTOM | DIR_TOP | DIR_FRONT | DIR_BACK
};
typedef size_t DirBC;

template <typename Map, typename Matrix, typename MultiVector>
class Problem : public Teuchos::Describable {
 public:
  Problem(Teuchos::ParameterList& list)
    : list_(list) {
    SetBoundary();
  };
  Problem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
    : list_(list) {
    Map_ = map;
    SetBoundary();
  };
  virtual ~Problem() {}

  virtual Teuchos::RCP<Matrix> BuildMatrix() = 0;
  virtual Teuchos::RCP<MultiVector> BuildCoords() {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Coordinates construction is not implemented for this problem");
  }
  virtual Teuchos::RCP<MultiVector> BuildNullspace() {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Nullspace construction is not implemented for this problem");
  }

  // Get methods
  Teuchos::RCP<const Map> getMap() const { return Map_; }
  Teuchos::RCP<const Matrix> getMatrix() const { return A_; }
  Teuchos::RCP<const MultiVector> getNullspace() const { return Nullspace_; }
  Teuchos::RCP<const MultiVector> getCoords() const { return Coords_; }

  // Set methods
  Teuchos::RCP<const Map> setMap(const Teuchos::RCP<const Map>& map) { Map_ = map; }

 protected:
  Teuchos::ParameterList& list_;
  Teuchos::RCP<const Map> Map_;
  Teuchos::RCP<Matrix> A_;
  Teuchos::RCP<MultiVector> Nullspace_;
  Teuchos::RCP<MultiVector> Coords_;

  DirBC DirichletBC_;

 private:
  void SetBoundary() {
    DirichletBC_ = DIR_ALL;
    if (this->list_.get("left boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_LEFT;
    if (this->list_.get("right boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_RIGHT;
    if (this->list_.get("bottom boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_BOTTOM;
    if (this->list_.get("top boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_TOP;
    if (this->list_.get("front boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_FRONT;
    if (this->list_.get("back boundary", "Dirichlet") == "Neumann") DirichletBC_ ^= DIR_BACK;
  }
};

}  // namespace Xpetra

}  // namespace ADR

#endif  // ADR_PROBLEM_HPP
