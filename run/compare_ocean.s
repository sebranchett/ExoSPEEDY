#!/bin/sh
set -e

# check a sample of values in .txt output files
diff -w checked_ocean_output.txt test_ocean_output.txt
