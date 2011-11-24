#!/bin/bash

getFileHash() {
	hash=$(md5sum $1  | awk '{print $1}');
	echo $hash;
} 



CodedDirToVerify=../Coded
CodedDirFinal=coded_final

if [ "$1" = "small" ]; then
	#echo "small"
	Files_L="testPatterns_1024"
elif [ "$1" = "medium" ]; then
	#echo "medium"
	Files_L="ILikeEI"
else 
	#echo "no params";
	Files_L="Informatica_e_no_IPL---cartaz
	ILikeEI
	testPatterns_1024"
fi



for File in $Files_L
do

	FileFinal=$CodedDirFinal/${File}_dic_1024_4x4.coded
	FileToVerify=$CodedDirToVerify/${File}_dic_1024_4x4.coded
	
	if [ $(getFileHash $FileFinal) = $(getFileHash $FileToVerify) ]; then
		echo "MATCHES ===== ${File}_dic_1024_4x4 !";
	else 
		echo "DOES NOT MATCH ===== ${File}_dic_1024_4x4 !";
		exit 0;
	fi
	
done
