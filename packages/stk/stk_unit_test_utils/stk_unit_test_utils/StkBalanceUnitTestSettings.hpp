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

} }

#endif // STKBALANCEUNITTESTSETTINGS_HPP
