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
keyFingerprint="2048 b4:9c:c1:ef:95:12:df:54:a3:98:13:9e:db:68:fb:3f csiefer@sandia.gov (RSA)"
# socket query tool
socketCommand="/usr/sbin/ss"

shell = os.environ["SHELL"]
if shell == "/bin/bash":
   envCmd="export SSH_AUTH_SOCK="
elif shell == "/bin/tcsh":
   envCmd="printenv SSH_AUTH_SOCK "
else:
  print "Only bash and tcsh are supported."
  quit()

# Your username.
myName = pwd.getpwuid(os.getuid())[0]
(status,charlist)=commands.getstatusoutput(socketCommand + " -xl | grep -o '/tmp/ssh-[[:alnum:]]*/agent.[[:digit:]]*'")
# convert raw characters into list
agentList=[s.strip() for s in charlist.splitlines()]
agentFound=0
keyFound=0
for agent in agentList:
  #see if this is your agent by checking ownership of root of lock directory
  #check only the root, because if it's not yours, you can't see down into it.
  pieces=agent.split("/")
  rootDir = "/" + pieces[1] + "/" + pieces[2]
  st = os.stat(rootDir)
  dirOwner = pwd.getpwuid(st.st_uid).pw_name
  if dirOwner == myName:
    agentFound=1
    #your ssh agent has been found
    sshAgentCmd="SSH_AUTH_SOCK=" + agent + " ssh-add -l"
    (status,result)=commands.getstatusoutput(sshAgentCmd)
    keyList=[s.strip() for s in result.splitlines()]
    #check whether this key's fingerprint matches the desired key's
    for key in keyList:
      if key == keyFingerprint:
        keyFound=1
        #print "export SSH_AUTH_SOCK=" + agent
        print envCmd + agent
        break

#if no key matches, just use the last agent found
if keyFound == 0 and agentFound == 1:
  #print "export SSH_AUTH_SOCK=" + agent
  print envCmd + agent

