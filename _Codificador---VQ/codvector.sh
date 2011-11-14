#!/bin/bash

#InputFName=testPatterns_1024
#InputFName=tecnologias_emergentes02
InputDir_PGM=../Examples---PGM
DictsDir=../Dicts
CodedDir=../Coded
DictFName=dic_1024_4x4

Files_L="ILikeEI"

for File in $Files_L
do
	InputFName=$File
	InputFName_PGM=$InputFName.pgm
	echo "================================================================="
	echo ">> $DictFName"
	CodedFName=${InputFName}_${DictFName}.coded
	time ./codvector $InputDir_PGM/$InputFName_PGM $DictsDir/$DictFName $CodedFName
	ls -l $CodedFName
	mv $CodedFName $CodedDir

done
