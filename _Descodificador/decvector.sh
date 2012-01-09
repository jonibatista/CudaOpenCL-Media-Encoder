nputDir_PGM=../Images
DictsDir=../Dicts
CodedDir=../Coded

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

DictFName=dic_1024_4x4
#DictFName=dic_1024_2x2

for File in $Files_L
do
	BaseFName=$File
	OriginalFName=${BaseFName}.pgm

	# 4x4
	./decvector $CodedDir/${BaseFName}_${DictFName}.coded $DictsDir/$DictFName ../Decoded/${BaseFName}.rec.pgm $InputDir_PGM/$BaseFName.pgm

done
