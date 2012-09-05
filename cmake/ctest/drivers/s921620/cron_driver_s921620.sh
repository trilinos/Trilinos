
export MAILTO="shylu-regression@software.sandia.gov"
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80/"
#export TDD_CTEST_TEST_TYPE=Nightly
export HTTP_PROXY="http://sonproxy.sandia.gov:80/"
/usr/local/bin/git pull
#env Trilinos_PACKAGES=ShyLU CTEST_TEST_TYPE=Nightly ../cron_driver.py
env Trilinos_PACKAGES=ShyLU ../cron_driver.py
