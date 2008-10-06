#ifndef INTREPID_CUBATURE_SPARSE_HPP
#define INTREPID_CUBATURE_SPARSE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_CubatureSparseHelper.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Teuchos_TestForException.hpp"
//#include <vector.h>

namespace Intrepid{

//Need to talk to Denis about
int INTREPID_MAX_CUBATURE_DEGREE_SPARSE2D = 59;
int INTREPID_MAX_CUBATURE_DEGREE_SPARSE3D = 57;

template<class Scalar, int dimension_>
class CubatureSparse : public Intrepid::Cubature<Scalar> {
  private:
  int level_;
  int numPoints_;
  const int degree_;
  
  public:

  ~CubatureSparse() {}


  CubatureSparse(const int degree);

  virtual void getCubature(int &                            numCubPoints,
                           Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const;

  virtual void getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const;

  virtual int getNumPoints() const;

  ECell getCellType() const;

  virtual int getAccuracy() const;

}; // end class CubatureSparse 


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureSparseDef.hpp>

#endif
