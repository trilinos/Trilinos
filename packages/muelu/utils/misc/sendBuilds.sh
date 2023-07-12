#!/bin/bash

# Usage: sendBuilds.sh FILENAME EMAIL
# Dependencies: mail

# check inputs
if [ $# -ne 2 ]; then
  echo "Invalid number of arguments given!"
  echo "Usage: sendBuilds.sh FILENAME EMAIL"
  echo "   FILENAME the attachment to send"
  echo "   EMAIL the email address to send to"
  echo "Exiting..."
  exit 1
fi

# silently check if mail is available
which mail > /dev/null
if [ $? -ne 0 ]; then
  echo "Program mail not found in PATH!"
  echo "Exiting..."
  exit 1
fi

# check if the file to attach is valid
if [ -f "$1" ]; then
  # send file as attachment to the specified address
  echo "Sending $1 to address $2..."
  echo "PR Results attached" | mail -a $1 -s "PR Results: $1" $2
else
  echo "Invalid file specified!"
  echo "File $1 does not exist!"
  echo "Exiting..."
  exit 1
fi

echo "Done!"
