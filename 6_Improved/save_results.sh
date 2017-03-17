#!/bin/bash
message=$1
if [[ -n "$message" ]]; then
    mkdir ../Tests/"$message"/
    mv U* ../Tests/"$message"/
    ls ../Tests/"$message"/
else
    mv U* ../Tests/
    ls ../Tests/
fi
git add *
git commit -m "return results"
git push origin master
ls
