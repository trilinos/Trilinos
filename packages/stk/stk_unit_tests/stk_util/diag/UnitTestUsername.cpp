#include <gtest/gtest.h>
#include <stk_util/diag/Platform.hpp>

#include <string>
#include <cstdlib>

TEST(StkEnv, getenv)
{
  char* env_user_value = new char[64];
  sprintf(env_user_value, "USER=");
  putenv(env_user_value);

  char*env_user = std::getenv("USER");

  //getenv can return non-null but empty strings...
  if (env_user) {
    std::string str_env_user(env_user);
    if (str_env_user.empty()) {
      env_user = NULL;
    }
  }

  char*expected_env_user = NULL;

  EXPECT_EQ(expected_env_user, env_user);

  sprintf(env_user_value, "USER=the_man");
  putenv(env_user_value);

  env_user = std::getenv("USER");
  std::string str_env_user(env_user);

  std::string expected_user("the_man");

  bool they_match = expected_user == str_env_user;
  EXPECT_TRUE(they_match);

  delete [] env_user_value;
}

