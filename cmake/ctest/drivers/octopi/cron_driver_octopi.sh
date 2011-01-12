
export MAILTO="zoltan-regression@software.sandia.gov"
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80/"
export TDD_CTEST_TEST_TYPE=Nightly
export HTTP_PROXY="http://sonproxy.sandia.gov:80/"
export PATH=/opt/lam714-gcc346-pure/bin:/usr/local/bin:$PATH
/opt/lam714-gcc346-pure/bin/lamboot
env Trilinos_PACKAGES=Zoltan CTEST_TEST_TYPE=Nightly ../cron_driver.py
