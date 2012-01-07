#!/bin/bash

#InputFName=testPatterns_1024
#InputFName=tecnologias_emergentes02
InputDir_PGM=../Images
DictsDir=../Dicts
CodedDir=../Coded
DictFName=dic_1024_4x4


if [ "$1" = "small" ]; then
	#echo "small"
	Files_L="testPatterns_1024"
elif [ "$1" = "medium" ]; then
	#echo "medium"
	Files_L="ILikeEI"
elif [ "$1" = "large" ]; then
	#echo "large"
	Files_L="Informatica_e_no_IPL---cartaz"
else 
	#echo "no params";
	Files_L="Informatica_e_no_IPL---cartaz
	ILikeEI
	testPatterns_1024"
fi


for File in $Files_L
do
	InputFName=$File
	InputFName_PGM=$InputFName.pgm
	echo "================================================================="
	echo ">> $DictFName"
	CodedFName=${InputFName}_${DictFName}.coded
	time ./codvector -i $InputDir_PGM/$InputFName_PGM -d $DictsDir/$DictFName -o $CodedFName
	ls -l $CodedFName
	mv $CodedFName $CodedDir

done
