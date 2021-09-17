#!/bin/bash

TASK=Sandbox
TEST=test/${TASK}
SRC1=ExtremaFinding
SRC2=numerical_recipes
SRC3=ObjectiveFunctions
SRC4=SimplicityAccuracy

g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} ${TEST}/test.cpp -Isrc/${SRC1} src/${SRC1}/*.cpp -Isrc/${SRC2} src/${SRC2}/*.cpp -Isrc/${SRC3} src/${SRC3}/*.cpp -Isrc/${SRC4} src/${SRC4}/*.cpp 

./bin/${TASK}