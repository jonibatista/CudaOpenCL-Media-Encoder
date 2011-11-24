#!/bin/bash

getFileHash() {
	hash=$(md5sum $1  | awk '{print $1}');
	echo $hash;
} 



CodedDirToVerify=../Coded
CodedDirFinal=coded_final


Files_L="Informatica_e_no_IPL---cartaz
ILikeEI
testPatterns_1024"


for File in $Files_L
do

	FileFinal=$CodedDirFinal/${File}_dic_1024_4x4.coded
	FileToVerify=$CodedDirToVerify/${File}_dic_1024_4x4.coded
	
	#echo $(getFileHash $File);
	
	if [ $(getFileHash $FileFinal) = $(getFileHash $FileToVerify) ]; then
		echo "MATCHES ===== ${File}_dic_1024_4x4 !";
	else 
		echo "DOES NOT MATCH ===== ${File}_dic_1024_4x4 !";
	fi
	
	
	
done
