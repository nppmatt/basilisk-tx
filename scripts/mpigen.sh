#!/usr/bin/env bash

# Possible bug in Basilisk compiler, where the -source flag does not consider directory structure.
# Workaround is to just have a shell cd into the src directory and invoke the translation there.
cd src
qcc -source -D_MPI=1 $1

