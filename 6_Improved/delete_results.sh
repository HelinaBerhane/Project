#!/bin/bash
message=$1
if [[ -n "$message" ]]; then
    rm *"$message".txt
else
    rm U*
fi
ls
