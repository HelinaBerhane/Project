#!/bin/bash
message=$1
git add *
if [[ -n "$message" ]]; then
    git commit -m "$message"
else
    git commit -m "testing"
fi
git push origin master
#
