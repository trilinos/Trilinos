
export MAILTO="zoltan-regression@software.sandia.gov"
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80/"
export HTTP_PROXY="http://sonproxy.sandia.gov:80/"
export PATH=/opt/lam714-gcc346-pure/bin:/usr/local/bin:$PATH
env Trilinos_PACKAGES=Zoltan CTEST_TEST_TYPE=Experimental ../cron_driver.py
