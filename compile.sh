#!/bin/bash

TASK=ExtremaFinding

g++ -std=c++11 -Wall -O2 -static -o ${TASK} grader.cpp ${TASK}.cpp
