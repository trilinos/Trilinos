#ifndef INTREPID_CUBATURE_GEN_SPARSE_HPP
#define INTREPID_CUBATURE_GEN_SPARSE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_CubatureSparseHelper.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid{

//Need to talk to Denis about
 const int INTREPID_MAX_CUBATURE_DEGREE_GENSPARSE = 17;


template<class Scalar, int dimension_>
class CubatureGenSparse : public Intrepid::Cubature<Scalar> {
  private:
  int numPoints_;
  const int degree_;
  SGNodes<Scalar, dimension_> grid;
  
  public:

  ~CubatureGenSparse() {}


  CubatureGenSparse(const int degree);

  virtual void getCubature(int &                            numCubPoints,
                           Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const;

  virtual void getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const;

  virtual int getNumPoints() const{return numPoints_;};

  ECell getCellType() const{return CELL_HEX;};

  virtual int getAccuracy() const{return 0;};

}; // end class CubatureSparse 


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureGenSparseDef.hpp>

#endif
