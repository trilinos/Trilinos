#!/bin/bash

# Run top in batch mode (-b) to update every 15 minutes ('-b 900'
# seconds) repreated -n 98 times so it will cover an entire day
# (i.e. 24*60/15=100)

top -b -d 900 -n 98
