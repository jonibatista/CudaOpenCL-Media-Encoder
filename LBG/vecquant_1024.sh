#!/bin/bash

Size=1024

DirPGM=../Training_pgm/pgm_1024

DicFName=dic_${Size}_4x4
./vecquant $Size 4 4 $DicFName 1  $DirPGM/t_airport_1024.pgm $DirPGM/t_pentagon_1024.pgm $DirPGM/t_tileRoof_1024.pgm $DirPGM/man_1024.pgm
##ls -l $DicFName

DicFName=dic_${Size}_2x2
./vecquant $Size 2 2 $DicFName 1  $DirPGM/t_airport_1024.pgm $DirPGM/t_pentagon_1024.pgm $DirPGM/t_tileRoof_1024.pgm $DirPGM/man_1024.pgm

# Commented out: no execution time difference in coding between:
# - dic_512_4x4 
# - dic_512_2x2
# 
#
# DicFName=dic_512_2x2
# ./vecquant 512 2 2 $DicFName 1 ../PGM/Baboon.pgm ../PGM/flowers.pgm ../PGM/lena512.pgm ../PGM/peppers_gray.pgm
