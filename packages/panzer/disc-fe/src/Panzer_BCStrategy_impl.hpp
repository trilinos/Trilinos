// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_IMPL_HPP
#define PANZER_BCSTRATEGY_IMPL_HPP

template <typename EvalT>
panzer::BCStrategy<EvalT>::
BCStrategy(const panzer::BC& bc) :
  m_bc(bc) 
{ }

template <typename EvalT>
panzer::BCStrategy<EvalT>::~BCStrategy()
{ }

#endif
