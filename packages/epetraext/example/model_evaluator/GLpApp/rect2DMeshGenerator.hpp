
#include <iosfwd>

class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;
class Epetra_IntSerialDenseMatrix;

namespace GLpApp {

/** \brief Generate a simple rectangular 2D triangular mesh that is only
 * partitioned between processors in the <tt>y</tt> direction.
 *
 * ToDo: Finish documentation!
 */
void rect2DMeshGenerator(
  const int                      numProc
  ,const int                     procRank
  ,const double                  len_x
  ,const double                  len_y
  ,const int                     local_nx
  ,const int                     local_ny
  ,const int                     bndy_marker
  ,Epetra_IntSerialDenseVector   *ipindx_out
  ,Epetra_SerialDenseMatrix      *ipcoords_out
  ,Epetra_IntSerialDenseVector   *pindx_out
  ,Epetra_SerialDenseMatrix      *pcoords_out
  ,Epetra_IntSerialDenseMatrix   *t_out
  ,Epetra_IntSerialDenseMatrix   *e_out
  ,std::ostream                  *out
  ,const bool                    dumpAll
  );

} // namespace GLpApp
