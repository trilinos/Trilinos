#ifndef _MiniEM_Utils_hpp_
#define _MiniEM_Utils_hpp_

#include "Teuchos_RCP.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Panzer_NodeType.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Teko_Utilities.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"

namespace mini_em {

  void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op);

  void describeMatrix(const std::string & s,const Thyra::LinearOpBase<double> & op,Teuchos::RCP<Teuchos::FancyOStream> out);

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node=Tpetra::Map<>::node_type>
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > get_Tpetra_CrsMatrix(const Thyra::LinearOpBase<double> & op);

  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::LinearOpBase<double> & op);

  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::DiagonalLinearOpBase<double> & op, const Epetra_Comm& comm);

}
#endif
