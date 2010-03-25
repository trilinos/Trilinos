#!/bin/sh

# This commit backs out a bunch of bad commits that were made to
# Zoltan on 2010/03/13.  These can be reverted again with a command.

# A) Reverse the commits

eg cherry-pick -R 4e851528cbbe24e81c87de3a7f7763761cb558ee
eg cherry-pick -R 81c6073d949b6336034bbd358ca00a413f972bcb
eg cherry-pick -R 96014552a73d5083a1bcfb4b7442d00045f779fe
eg cherry-pick -R e2219ac5a7c510d1a66bba39a0a016242bfe7c16
eg cherry-pick -R 8b686a930bc0649f610c2d8c4c9724574986816f
eg cherry-pick -R 20ebe202a40af04d28f2d7e321b4d30ef927f453
eg cherry-pick -R 1f8c8ab86b811200f097615c4c9d6bc081b2cca2
eg cherry-pick -R 1f8c8ab
eg cherry-pick -R -m 1 4b9442e
# NOTE: Above, you have to specify the parent and it is 1 (the master
# branch) in this case

# B) Squash the commits into one

eg squash

# Where you get here you should be ready to test the commits that have
# been reverted.
