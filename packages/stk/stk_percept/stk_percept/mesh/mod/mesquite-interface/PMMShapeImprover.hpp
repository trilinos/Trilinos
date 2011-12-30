/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMShapeImprover_hpp
#define PMMShapeImprover_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <InstructionQueue.hpp>
#include <Settings.hpp>
#include <ShapeImprovementWrapper.hpp>

#include <IdealWeightInverseMeanRatio.hpp> 
#include <QualityAssessor.hpp>
#include <TerminationCriterion.hpp>

#include <TQualityMetric.hpp>
#include <TShapeB1.hpp>
#include <IdealShapeTarget.hpp>

#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>


namespace stk {
  namespace percept {

    using namespace Mesquite;

    class PMMShapeImprover
    {
    public:
      PMMShapeImprover() {}

      static int count_invalid_elements(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain)
      {
        MsqError err;
        InstructionQueue q;
  
        // Set up barrier metric to see if mesh contains inverted elements
        TShapeB1 mu_b;
        IdealShapeTarget w_ideal;
        TQualityMetric barrier( &w_ideal, &mu_b );
  
        // Check for inverted elements in the mesh
        QualityAssessor inv_check( &barrier );
        //inv_check.disable_printing_results();
        q.add_quality_assessor( &inv_check, err );  MSQ_ERRZERO(err);
        Settings settings;
        q.run_common( &mesh, 0, &domain, &settings, err ); MSQ_ERRZERO(err);
        //q.remove_quality_assessor( 0, err ); MSQ_ERRZERO(err);
        const QualityAssessor::Assessor* inv_b = inv_check.get_results( &barrier );
        int num_invalid = inv_b->get_invalid_element_count();
        return num_invalid;
      }

      void run(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain, bool always_smooth=true, bool debug=false)
      {
        if (debug)
          {
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            Mesquite::MsqDebug::enable(3);
          }
        Mesquite::MsqError mErr;
        int num_invalid = 0;
        bool check_quality=true;
        if (check_quality)
          {
            num_invalid = count_invalid_elements(mesh, domain);
            std::cout << "tmp srk PMMShapeImprover num_invalid before= " << num_invalid << std::endl;
          }

        if (num_invalid || always_smooth)
          {
            Mesquite::ShapeImprovementWrapper siw(mErr);
            //siw.set_iteration_limit(1);

            siw.run_instructions(&mesh, &domain, mErr);

            if (check_quality)
              {
                num_invalid = count_invalid_elements(mesh, domain);
                std::cout << "tmp srk PMMShapeImprover num_invalid after= " << num_invalid << std::endl;
              }
          }
      }
    };

  }
}

#endif
#endif
#endif
