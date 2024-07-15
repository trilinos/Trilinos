// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOMAIN_INTERFACE_HPP
#define PANZER_EVALUATOR_DOMAIN_INTERFACE_HPP

namespace panzer {

  class Workset;

  /** \brief Mix-in interface to support cell "domains" in panzer.
   * 
   *  This class adds support for cell domains into evaluators. This is commonly used in DG discretizations.
   */
  class DomainEvaluator {

  public:

    /// Domain types supported by worksets. 
    enum DomainType : int {
        OWNED=0,    /// All Owned cells for the workset on the MPI process
        GHOST=1,    /// All Ghosted cells for the workset on the MPI process
        REAL=2,     /// All Owned and Ghosted cells for the workset on the MPI process
        VIRTUAL=3,  /// All virtual cells for the workset on the MPI process
        EXTERNAL=4, /// All ghost and virtual cells for the workset on the MPI process
        ALL=5       /// All OWNED, GHOSTED and VIRTUAL cells for the workset on the MPI process
        };

    /**
     * \brief Constructor
     *
     * \param[in] domain (optional) Cell domain to iterate over (defaults to ALL)
     */
    DomainEvaluator(DomainType domain=ALL);

    /**
     * \brief Default destructor
     */
    virtual ~DomainEvaluator() = default;

    /**
     * \brief Set the domain for the evaluator
     *
     * \param[in] domain Domain to set
     */
    void setDomain(const DomainType domain);
    
    /**
     * \brief Get the domain for the evaluator
     *
     */
    DomainType getDomain();

    /**
     * \brief Returns the starting cell for the specified domain for a given workset
     *
     * Note: the loop would look like:
     *
     * for(int cell = cell_start_index(workset); cell < cell_end_index(workset); ++cell){do something...}
     *
     * \param[in] workset Workset describing data layout
     *
     * \return Starting index for cell domain
     */
    virtual int cellStartIndex(const panzer::Workset & workset) const;

    /**
     * \brief Returns the non-inclusive end cell for the specified domain for a given workset
     *
     * Note: the loop would look like:
     *
     * for(int cell = cell_start_index(workset); cell < cell_end_index(workset); ++cell){do something...}
     *
     * \param[in] workset Workset describing data layout
     *
     * \return End index for cell domain
     */
    virtual int cellEndIndex(const panzer::Workset & workset) const;

  private:

    /// Domain for this evaluator
    DomainType domain_;

  };

}

#endif
