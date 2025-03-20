#!/bin/bash

# Usage: ./test.sh [--regold] [test-name]
# Example:
#   Run all tests:
#     ./test.sh
#   Regold all tests:
#     ./test.sh --regold
#   Run one test:
#     ./test.sh tests/gold-45-sas-h-nb.txt
#   Regold one test:
#     ./test.sh --regold tests/gold-45-sas-h-nb.txt
./test-align.sh "$@"
