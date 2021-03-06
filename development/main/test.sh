#!/bin/bash

TASK=main
TEST=development/${TASK}
SRC=ExtremaFinding

g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} ${TEST}/test.cpp -Isrc/numerical_recipes src/numerical_recipes/*.cpp -Isrc/${SRC} src/${SRC}/*.cpp

./bin/${TASK}
