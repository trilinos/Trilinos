#!/usr/bin/env python

import commands
import os
import re
import sys

if len(sys.argv) != 3:
  raise SystemExit("Syntax:\n  %s oldrev newrev" % sys.argv[0])

oldrev=sys.argv[1]
newrev=sys.argv[2]

output=commands.getoutput("git diff --name-only %s %s" % (oldrev, newrev))
dirschanged = [os.path.dirname(filename) for filename in output.splitlines()]
defaultEmail = commands.getoutput("eg config --get hooks.mailinglist").strip()
#print "defaultEmail =", defaultEmail

dirs = {}.fromkeys(dirschanged, 1)
#print "dirs =", dirs

# Create a list of (regex, email) in the order they are in the file.  That
# way, the first regex that matches will be used to select the email address.

f = open(os.path.join(os.path.dirname(sys.argv[0]),'dirs_to_emails'), 'r')

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
for dir in dirs:
  #print "\ndir =", dir
  found = False
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
