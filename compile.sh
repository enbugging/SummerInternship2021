#!/bin/bash

TASK=ExtremaFinding
TASK1=GlobalMinimumFinder
TASK2=LocalMinimaFinder

g++ -c -I ./src/numerical_recipes ./src/numerical_recipes/*.cpp
g++ -c -I ./src/ExtremaFinding ./src/ExtremaFinding/*.cpp
g++ -std=c++11 -Wall -O2 -static -o grader.cpp -I .
#rm *.o
#g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} grader.cpp -Ibin/numerical_recipes.o -I${TASK} # ${TASK}/${TASK1}.cpp ${TASK}/${TASK2}.cpp 