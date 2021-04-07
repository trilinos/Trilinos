#!/usr/bin/env python
import os
import pwd
import commands
import re

#
# Utility to find out if you have an ssh-agent running that is holding your
# private key.  To use this in bash:
#
#     eval $(python ./setSetSshEnv.python)
#
# It assumes that ssh creates files of the form /tmp/ssh-Abcdefg12345/agent.12345 .
#
# TODO 1: It would be better to skip the socket query and instead look directly for the ssh-agent lock files.

# Fingerprint of identity that you will use. You can find it with "ssh-add -l".
keyFingerprint = "4096 db:7b:ce:13:b4:88:64:59:b2:e9:b8:17:f0:a6:97:d5 id_rsa_geminga (RSA)"
# socket query tool
socketCommand="/usr/sbin/ss"

shell = os.environ["SHELL"]
if shell == "/bin/bash" or shell == "/bin/sh":
   envCmd="export SSH_AUTH_SOCK="
elif shell == "/bin/tcsh" or shell == "/bin/csh":
   envCmd="setenv SSH_AUTH_SOCK "
else:
  print "Only bash, csh, and tcsh are supported."
  quit()

# Your username.
userid = pwd.getpwuid(os.getuid())[0]
[status,charlist]=commands.getstatusoutput(socketCommand + " -xl | grep -o '/tmp/ssh-[[:alnum:]]*/agent.[[:digit:]]*'")
# convert raw characters into list
agentList=[s.strip() for s in charlist.splitlines()]
agentFound=0
keyFound=0
myagent = ""
for agent in agentList:
  # See if this is your agent by checking ownership of root of lock directory
  # Check only the root, because if it's not yours, you can't see down into it.
  pieces=agent.split("/")
  rootDir = "/" + pieces[1] + "/" + pieces[2]
  # JJH: On redsky, the socket command returned nonexistent directories
  #      So I check for existence first to avoid an exception when calling os.stat
  #      on a nonexistent directory.
  if os.path.isdir(rootDir):
    st = os.stat(rootDir)
    dirOwner = pwd.getpwuid(st.st_uid).pw_name
    if dirOwner == userid:
      agentFound=1
      myagent=agent
      # Your ssh agent has been found
      sshAgentCmd="SSH_AUTH_SOCK=" + myagent + " ssh-add -l"
      [status,result]=commands.getstatusoutput(sshAgentCmd)
      keyList=[s.strip() for s in result.splitlines()]

      # Check whether this key's fingerprint matches the desired key's
      for key in keyList:
        if keyFingerprint in key:
          keyFound=1
          print envCmd + myagent
          break

# If no key matches, just use the last agent found
if keyFound == 0 and agentFound == 1:
  #print "export SSH_AUTH_SOCK=" + agent
  print envCmd + myagent
