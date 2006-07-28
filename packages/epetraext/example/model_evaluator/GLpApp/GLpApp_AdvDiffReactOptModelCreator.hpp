#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

namespace GLpApp {

/** \brief A utility class for creating an <tt>AdvDiffReactOptModelCreator/tt>
 * object by reading from the command-line.
 */
class AdvDiffReactOptModelCreator {
public:

  /** \brief . */
  AdvDiffReactOptModelCreator();

  /** \brief . */
  void setupCLP( Teuchos::CommandLineProcessor *clp );

  /** \brief . */
  Teuchos::RefCountPtr<AdvDiffReactOptModel>
  createModel(
    const Teuchos::RefCountPtr<const Epetra_Comm>     &comm
    ,std::ostream                                     *out  = NULL
    ) const;

private:

    double              len_x_;
    double              len_y_;
    int                 local_nx_;
    int                 local_ny_;
    std::string         geomFileBase_;
    int                 np_;
    bool                normalizeBasis_;
    double              reactionRate_;
    double              beta_;
    double              x0_;
    double              p0_;
  
};

} // namespace GLpApp
