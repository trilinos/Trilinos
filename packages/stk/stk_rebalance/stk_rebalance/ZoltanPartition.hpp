/*----------------------------------------------------------------------*/
/*                                                                      */
/*       author: Jonathan Scott Rath                                    */
/*      author2: Michael W. Glass   (DEC/2000)                          */
/*     filename: ZoltanPartition.h                                      */
/*      purpose: header file for stk toolkit zoltan methods             */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*    Copyright 2001,2010 Sandia Corporation.                           */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a         */
/*    non-exclusive license for use of this work by or on behalf        */
/*    of the U.S. Government.  Export of this program may require       */
/*    a license from the United States Government.                      */
/*----------------------------------------------------------------------*/

// Copyright 2001 Sandia Corporation, Albuquerque, NM.

#ifndef stk_rebalance_ZoltanPartition_hpp
#define stk_rebalance_ZoltanPartition_hpp

#include <utility>
#include <vector>
#include <string>

#include <Teuchos_ParameterList.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_rebalance/GeomDecomp.hpp>

//Forward declaration for pointer to a Zoltan structrue.
struct Zoltan_Struct;

namespace stk {
namespace rebalance {

class Zoltan : public GeomDecomp {

public:

  typedef std::map<std::string,std::string> Parameters;

  static const std::string zoltan_parameters_name();
  static const std::string default_parameters_name();

  static void init_default_parameters();


  /**
   * Constructor
   */

  static Zoltan create_default(ParallelMachine pm,  const unsigned ndim, const Teuchos::ParameterList & rebal_region_parameters);

  explicit Zoltan(ParallelMachine pm, const unsigned ndim, const std::string &Parameters_Name=default_parameters_name());

  void init(const std::vector< std::pair<std::string, std::string> >
            &dynamicLoadRebalancingParameters);

  static double init_zoltan_library();
  /**
   * Destructor
   */

  virtual ~Zoltan();

  /** Name Conversion Functions.
   * Long friendly string prarameters
   * need to be converted into short Zoltan names.  These
   * two functions do that.  Merge_Default_Values should be called
   * first because it merges with the long names.  Convert_Names_and_Values
   * should then be called second to convert long names to short names.
   */
  static void merge_default_values   (const Teuchos::ParameterList &from,
                                      Teuchos::ParameterList &to);
  static void convert_names_and_values (const Teuchos::ParameterList &from,
                                        Teuchos::ParameterList &to);

  /**
   * Register SIERRA Framework Zoltan call-back functions */
  int register_callbacks();

  virtual int determine_new_partition (bool & RebalancingNeeded);

  virtual int get_new_partition(std::vector<mesh::EntityProc> &new_partition)
  { throw std::runtime_error("stk::rebalance::Zoltan::get_new_partition not yet implemented."); }

  /**
   * Evaluate the performance/quality of dynamic load rebalancing
   */
  int evaluate ( int   print_stats,
                 int   *nobj,
                 double  *obj_wgt,
                 int   *ncuts,
                 double  *cut_wgt,
                 int   *nboundary,
                 int   *nadj         );

  /**
   * Decomposition Augmentation
   */
  virtual int point_assign( double    *coords,
                            int  *proc ) const;

  virtual int box_assign ( double min[],
                           double max[],
                           std::vector<int> &procs) const;

  /**
   * Inline functions to access private data
   */

  double zoltan_version()  const;
  const std::string & parameter_entry_name() const;

  Zoltan_Struct * zoltan() {
    return zoltan_id;
  }
  const Zoltan_Struct * zoltan() const {
    return zoltan_id;
  }
  unsigned spatial_dimension() const {
    return m_spatial_dimension_;
  }
private:
  /** Zoltan load balancing struct       */
  struct    Zoltan_Struct *zoltan_id;

  const     unsigned       m_spatial_dimension_;
  /** Name that was used to initialize this Zoltan_Struct
   * if the parameter constructor was used.
   */
  std::string  parameter_entry_Name;

  static const std::string zoltanparametersname;
  static const std::string defaultparametersname;
  Parameters               m_default_parameters_;


};

}
} // namespace stk

#endif
