#!/bin/bash

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

function check() {
  local error_code
  local this_test="$1"
  local gold="tests/${PREFIX}.${this_test}.gold"
  shift

  if [ "$TEST_NAME" != "all" -a "$TEST_NAME" != "${this_test}" ]; then
    return
  fi

  echo "    check_gold $gold matches $@"
  local tmp="/tmp/icdiff.output"
  $INVOCATION "$@" --columns=80 &> "$tmp"
  error_code="$?"

  if $REGOLD; then
    if [ -e "$gold" ] && diff "$tmp" "$gold" > /dev/null; then
      echo "Did not need to regold $gold"
    else
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
