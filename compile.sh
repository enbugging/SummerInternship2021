#!/bin/bash

TASK=ExtremaFinding
TASK1=GlobalMinimumFinder
TASK2=LocalMinimaFinder

g++ -std=c++11 -Wall -O2 -static -o bin/${TASK} grader.cpp -Isrc/numerical_recipes src/numerical_recipes/*.cpp -Isrc/${TASK} src/${TASK}/*.cpp
