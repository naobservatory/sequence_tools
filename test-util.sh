#!/bin/bash

# Test framework.
#
# Each test should look like either:
#
#    check ...
#
# Or:
#
#    check_pipeline ...
#
# In either case, the framework will verify that running the command produces
# the expected "gold" output.  To update the gold outputs to match current
# behavior, run with --regold.
#
# By default runs all the tests, or you can specify the name of a specific test
# to run.

if [ "$#" -gt 0 -a "$1" = "--regold" ]; then
  REGOLD="true"
  shift
else
  REGOLD="false"
fi

TEST_NAME=all
if [ "$#" -gt 0 ]; then
  TEST_NAME="$1"
  shift
fi

if [ "$#" != 0 ]; then
  echo "Usage: '$0 [--regold] [test-name]'"
  exit 1
fi

function fail() {
  echo "FAIL"
  exit 1
}

# Usage: check test-name arg1 arg2 arg3
#
# Verify that the command we're testing (stored in global variable INVOCATION
# along with standard arguments) produces the expected output when run on the
# arguments.
function check() {
  local this_test="$1"
  shift
  check_pipeline "$this_test" "$INVOCATION $@"
}

# Usage: check-pipeline test-name 'foo | bar | baz'
#
# Verify that the pipeline we're testing (argument 2) produces the expected output.
function check_pipeline() {
  local error_code
  local this_test="$1"
  local gold="tests/${PREFIX}.${this_test}.gold"
  shift

  if [ "$TEST_NAME" != "all" -a "$TEST_NAME" != "${this_test}" ]; then
    return
  fi

  echo "    check_gold $gold matches $@"
  local tmp="/tmp/icdiff.output"
  bash -c "$*" &> "$tmp"
  error_code="$?"

  if $REGOLD; then
    if [ -e "$gold" ] && diff "$tmp" "$gold" > /dev/null; then
      echo "Did not need to regold $gold"
    else
      if [ -e "$gold" ]; then
        echo "Previous gold:"
        echo "--------------"
        cat "$gold"
        echo "Proposed gold:"
        echo "-----_--------"
      fi
      cat "$tmp"
      echo
      read -p "Is this correct? y/n > " -n 1 -r
      echo
      if [[ "$REPLY" =~ ^[Yy]$ ]]; then
        mv "$tmp" "$gold"
        echo "Regolded $gold."
      else
        echo "Did not regold $gold."
      fi
    fi
    return
  fi

  if ! diff "$gold" "$tmp"; then
    echo "Got: ($tmp)"
    cat "$tmp"
    echo "Expected: ($gold)"
    cat "$gold"
    fail
  fi
}
