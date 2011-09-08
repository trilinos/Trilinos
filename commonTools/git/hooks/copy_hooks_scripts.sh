#!/bin/sh

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/copy_hooks_scripts.sh/\1/g"`

#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

cp $_SCRIPT_DIR/post-receive .
cp $_SCRIPT_DIR/post-receive-email .
cp $_SCRIPT_DIR/get_recipients.py .
cp $_SCRIPT_DIR/update_push_log.py .
