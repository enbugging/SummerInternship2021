#!/bin/bash

TASK=ExtremaFinding
TASK1=GlobalMinimumFinder
TASK2=LocalMinimaFinder

g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} grader.cpp ${TASK}/${TASK1}.cpp #${TASK}/${TASK2}.cpp
./bin/${TASK}
