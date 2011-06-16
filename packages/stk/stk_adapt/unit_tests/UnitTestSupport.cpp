/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <unit_tests/UnitTestSupport.hpp>

namespace stk {
  namespace adapt {
    namespace unit_tests {

      /// CONFIGURATIONS
      ///  you can choose to use regression testing or not (see always_do_regression_tests)


      void UnitTestSupport::save_or_diff(PerceptMesh& eMesh, std::string filename, int option )
      {
        if (always_do_regression_tests || option == 1)
          {
            unsigned p_size = eMesh.getParallelSize();
            unsigned p_rank = eMesh.getParallelRank();
            std::string par_filename = filename;
            if (p_size > 1)
              {
                par_filename = filename+"."+toString(p_size)+"."+toString(p_rank);
              }
            if (Util::file_exists(par_filename))
              {
                int spatialDim = eMesh.getSpatialDim();

                PerceptMesh eMesh1(spatialDim);
                eMesh.saveAs("./tmp.e");
                eMesh1.openReadOnly("./tmp.e");

                PerceptMesh eMesh_gold(spatialDim);
                eMesh_gold.openReadOnly(filename);
                //eMesh_gold.printInfo("gold copy: "+filename, 2);
                //eMesh1.printInfo("compare to: "+filename, 2);
                {
                  std::string diff_msg = "gold file diff report: "+filename+" \n";
                  bool print_during_diff = false;
                  bool diff = PerceptMesh::mesh_difference(eMesh1, eMesh_gold, diff_msg, print_during_diff);
                  if (diff)
                    {
                      //std::cout << "tmp writing and reading to cleanup parts" << std::endl;

                      // write out and read back in to cleanup old parts
                      eMesh1.saveAs("./tmp.e");
                      PerceptMesh eMesh2(spatialDim);
                      eMesh2.openReadOnly("./tmp.e");
                      //std::cout << "tmp done writing and reading to cleanup parts" << std::endl;
                      bool diff_2 = PerceptMesh::mesh_difference(eMesh2, eMesh_gold, diff_msg, print_during_diff);
                      //std::cout << "tmp diff_2= " << diff_2 << std::endl;
                      diff = diff_2;
                      if (diff_2)
                        {
                          //bool diff_3 =
                            PerceptMesh::mesh_difference(eMesh2, eMesh_gold, diff_msg, true);
                        }
                    }

                  STKUNIT_EXPECT_TRUE(!diff);
                }
              }
            else
              {
                eMesh.saveAs(filename);
              }
          }
        else
          {
            eMesh.saveAs(filename);
          }
      }




      bool UnitTestSupport::always_do_regression_tests = true;

    }
  }
}

