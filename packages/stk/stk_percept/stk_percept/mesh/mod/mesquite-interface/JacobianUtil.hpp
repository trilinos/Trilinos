/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef JacobianUtil_hpp
#define JacobianUtil_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <AveragingQM.hpp>
#include <Vector3D.hpp>
#include <MsqMatrix.hpp>
#include <Exponent.hpp>

#include <stk_percept/PerceptMesh.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    class JacobianUtil : public Mesquite::AveragingQM
    {

    public:
      typedef double Vec3D[3];
      Vec3D mCoords[4];

      enum { NELEM_TYPES = 10, NNODES_MAX = 8 };

      double   m_detJ[NNODES_MAX]; 
      MsqMatrix<3,3> m_J[NNODES_MAX];
      //MsqMatrix<3,3> m_dJ[NNODES_MAX][NNODES_MAX][3];
      //typedef MsqMatrix<3,3> dJdXn_i[NNODES_MAX][3];
      //static dJdXn_i m_dJ[NELEM_TYPES][NNODES_MAX];
      //static bool m_dJ_init;
      MsqMatrix<3,3> m_dMetric_dA[NNODES_MAX];
      double m_grad[NNODES_MAX][NNODES_MAX][3];
      int m_num_nodes;
      bool m_scale_to_unit;

      JacobianUtil() : 
        AveragingQM( QualityMetric::LINEAR ),
        m_num_nodes(0),
        m_scale_to_unit(false)
      {
      }

      bool operator()(double& averageJ, PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field,
                      const CellTopologyData * topology_data_in = 0 );

      bool grad_metric_util( PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field,
                             const CellTopologyData * topology_data );

    private:
      bool jacobian_matrix_3D(double &detJ, MsqMatrix<3,3>& A, const Vec3D *x);
      bool jacobian_matrix_2D(double &detJ, MsqMatrix<3,3>& A, const Vec3D *x);

      void grad_util_2d(const MsqMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util(const MsqMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);

    };
  }
}

#endif
#endif
