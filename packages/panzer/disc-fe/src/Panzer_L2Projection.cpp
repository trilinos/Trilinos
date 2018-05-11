// @HEADER
// @HEADER

#include "PanzerDiscFE_config.hpp"
#include "Panzer_L2Projection.hpp"
#include "Panzer_L2Projection_impl.hpp"

template class panzer::L2Projection<int,panzer::Ordinal64>;

// template
// Teuchos::RCP<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
// panzer::buildMassMatrix<int,panzer::Ordinal64>(const panzer::BasisDescriptor& basis,
//                                                const panzer::IntegrationDescriptor& integrationDescriptor,
//                                                const Teuchos::RCP<Teuchos::MpiComm<int>>& comm,
//                                                const Teuchos::RCP<const panzer::ConnManager<int,panzer::Ordinal64>> connManager,
//                                                const std::vector<std::string> elementBlockNames,
//                                                const Teuchos::RCP<panzer::WorksetContainer>& worksetContainer);
