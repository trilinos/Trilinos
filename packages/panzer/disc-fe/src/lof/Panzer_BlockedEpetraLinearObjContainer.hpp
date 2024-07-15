// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_BlockedEpetraLinearObjContainer_hpp__
#define __Panzer_BlockedEpetraLinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Epetra_Map.h"

#include <unordered_map>

namespace panzer {

/** Linear object container for Block operators, this
  * always assumes the matrix is square.
  */
class BlockedEpetraLinearObjContainer : public LinearObjContainer
                                      , public ThyraObjContainer<double> {
public:
   typedef Thyra::VectorBase<double> VectorType;
   typedef Thyra::LinearOpBase<double> CrsMatrixType;

   //! Make sure row and column spaces match up 
   bool checkCompatibility() const;

   virtual void clear();

   //! Put a particular scalar in the matrix
   void initializeMatrix(double value);

   void setMapsForBlocks(const std::vector<Teuchos::RCP<const Epetra_Map> > & blockMaps)
   { blockMaps_ = blockMaps; }

   Teuchos::RCP<const Epetra_Map> getMapForBlock(std::size_t i) const
   { return blockMaps_[i]; }

   inline void set_x(const Teuchos::RCP<VectorType> & in) { set_x_th(in); }
   inline Teuchos::RCP<VectorType> get_x() const { return get_x_th(); }

   inline void set_dxdt(const Teuchos::RCP<VectorType> & in) { set_dxdt_th(in); }
   inline Teuchos::RCP<VectorType> get_dxdt() const { return get_dxdt_th(); }

   inline void set_f(const Teuchos::RCP<VectorType> & in) { set_f_th(in); }
   inline Teuchos::RCP<VectorType> get_f() const { return get_f_th(); }

   inline void set_A(const Teuchos::RCP<CrsMatrixType> & in) { set_A_th(in); }
   inline Teuchos::RCP<CrsMatrixType> get_A() const { return get_A_th(); }

   // Inherited from LinearObjContainer
   virtual void initialize();

   // Inherited from ThyraObjContainer
 
   void set_x_th(const Teuchos::RCP<VectorType> & in) { x = in; }
   Teuchos::RCP<VectorType> get_x_th() const { return x; }

   void set_dxdt_th(const Teuchos::RCP<VectorType> & in) { dxdt = in; }
   Teuchos::RCP<VectorType> get_dxdt_th() const { return dxdt; }

   void set_f_th(const Teuchos::RCP<VectorType> & in) { f = in; }
   Teuchos::RCP<VectorType> get_f_th() const { return f; }

   void set_A_th(const Teuchos::RCP<CrsMatrixType> & in) { A = in; }
   Teuchos::RCP<CrsMatrixType> get_A_th() const { return A; }

private:
   Teuchos::RCP<VectorType> x, dxdt, f;
   Teuchos::RCP<CrsMatrixType> A;

   std::vector<Teuchos::RCP<const Epetra_Map> > blockMaps_;
};

}

#endif
