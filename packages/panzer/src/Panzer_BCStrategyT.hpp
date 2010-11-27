#ifndef PANZER_BCSTRATEGY_T_HPP
#define PANZER_BCSTRATEGY_T_HPP

template <typename EvalT>
panzer::BCStrategy<EvalT>::
BCStrategy(const panzer::BC& bc) :
  m_bc(bc) 
{ }

#endif
