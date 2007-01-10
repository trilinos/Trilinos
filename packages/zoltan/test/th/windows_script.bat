@echo off

echo "Windows nightly batch file for Acro tests"
echo ""

C:
chdir C:\home\wehart\src\acro-th

echo "Launching the acro-th/test_daemon file... (please be patient)"

C:\cygwin\bin\tcsh -c /cygdrive/c/home/wehart/src/acro-th/cron_script > C:\home\wehart\src\acro-th\cron.out 2>&1
