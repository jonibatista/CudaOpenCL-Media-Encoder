#!/bin/bash

InputDir_PGM=../Examples---PGM
DictsDir=../Dicts
CodedDir=../Coded

Files_L="
Informatica_e_no_IPL---cartaz
ILikeEI
testPatterns_1024"


DictFName=dic_1024_4x4
#DictFName=dic_1024_2x2

for File in $Files_L
do
	BaseFName=$File
	OriginalFName=${BaseFName}.pgm

	# 4x4
	./decvector $CodedDir/${BaseFName}_${DictFName}.coded $DictsDir/$DictFName ../Decoded/${BaseFName}.rec.pgm $InputDir_PGM/$BaseFName.pgm

done
