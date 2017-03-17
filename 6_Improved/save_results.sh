#!/bin/bash
message=$1
if [[ -n "$message" ]]; then
    mkdir ../Tests/"$message"/
    mv U* ../Tests/"$message"/
else
    mv U* ../Tests/
fi
ls
