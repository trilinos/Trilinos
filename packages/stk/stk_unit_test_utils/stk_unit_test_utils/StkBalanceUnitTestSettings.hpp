#ifndef STKBALANCEUNITTESTSETTINGS_HPP
#define STKBALANCEUNITTESTSETTINGS_HPP

#include "stk_balance/balance.hpp"

namespace stk
{
namespace unit_test_util
{

class StkBalanceUnitTestSettings : public stk::balance::StkBalanceSettings
{
public:
  StkBalanceUnitTestSettings() : StkBalanceSettings() {}

  ~StkBalanceUnitTestSettings() = default;

  virtual bool getEdgesForParticlesUsingSearch() const override { return true; }
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
StkBalanceUnitTestSettings : public stk::unit_test_util::StkBalanceUnitTestSettings {};

} // namespace simple_fields

} }

#endif // STKBALANCEUNITTESTSETTINGS_HPP
