#!/usr/bin/env python
# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Test, and Integrate System
#                 Copyright (2011) Sandia Corporation
#
#
# Copyright (2011) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


import commands
import os
import re
import sys

if len(sys.argv) != 4:
  raise SystemExit("Syntax:\n  %s oldrev newrev refname" % sys.argv[0])

oldrev=sys.argv[1]
newrev=sys.argv[2]
refname=sys.argv[3]
index=refname.rfind("/")
shortrefname=refname[index+1:].lower()

output=commands.getoutput("git diff --name-only %s %s" % (oldrev, newrev))
dirschanged = [os.path.dirname(filename) for filename in output.splitlines()]
defaultEmail = commands.getoutput("git config --get hooks.mailinglist").strip()
#print "defaultEmail =", defaultEmail

dirs = {}.fromkeys(dirschanged, 1)
#print "dirs =", dirs

# Create a list of (regex, email) in the order they are in the file.  That
# way, the first regex that matches will be used to select the email address.

#print "sys.argv =", sys.argv
dirs_to_email_file = os.path.join(os.path.dirname(sys.argv[0]),'dirs_to_emails')
#print "dirs_to_email_file =", dirs_to_email_file
f = open(dirs_to_email_file, 'r')

emails = []

for raw_line in f:
  line = raw_line.strip()
  #print "\nline = '"+line+"'"
  if line.startswith('#') or line == "":
    continue
  (regex, email) = line.split()
  emails.append((regex,email))

f.close()

#for regex_email in emails:
#  print "%s: %s" % regex_email

recipients = {}

found = False

#if this is a branch for a package then only send the email to that package.
#otherwise send it to the list of packages that were modified.
#this is to cut down on the number of "extra" emails that get sent out when
#someone merges master onto a package branch and then pushes.
if not re.match("^master", shortrefname) and not re.match("^trilinos-release", shortrefname):
  for pair in emails:
    index = pair[1].find("-")
    packageRE = pair[1][:index]
    if re.match(packageRE, shortrefname):
      recipients[pair[1]] = 1
      found = True
      break
else:
  for dir in dirs:
    #print "\ndir =", dir
    for regex_email in emails:
      (regex, email) = regex_email
      #print "\nregex =", regex
      #print "email =", email
      if re.match(regex, dir):
        recipients[email] = 1
        found = True
        #print "FOUND IT!"
        break

if not found:
  recipients[defaultEmail] = 1
  #print "\nNOT FOUND: Use default email!"

print ",".join(sorted(recipients.keys()))
