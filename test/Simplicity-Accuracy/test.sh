#!/bin/bash

TASK-'Simplicity-Accuracy'
TEST=test/${TASK}
SRC=ExtremaFinding

g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} ${TEST}/test.cpp -Isrc/numerical_recipes src/numerical_recipes/*.cpp -Isrc/${SRC} src/${SRC}/*.cpp

./bin/${TASK}
