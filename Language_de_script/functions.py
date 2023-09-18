#!/usr/bin/env python3 
import re

# Définition des fonctions
def getfasta(file): 
	path_split = file.split("/")
	file_extension = path_split[len(path_split)-1].split(".")[1]

	if (file_extension!="fasta" and file_extension!="fa") :
		print("\n---------------------------TP1_languages_script.py------------------------------\nError[1]: Not a fasta file.\n")
		exit(1)

	while True :
		try :
			nameHandle = open(file,'r')
			break
		except FileNotFoundError :
			print("\n---------------------------TP1_languages_script.py------------------------------\nError[2]: Can't find file.\n")
			exit(2)
		except OSError :
			print("\n---------------------------TP1_languages_script.py------------------------------\nError[3]: Can't open/read the file.\n")
			exit(3)

	fastas = {} 
	for line in nameHandle: 
		if line[0]=='>': 
			header = line[1:-1]
			fastas[header]=''
		else: 
			fastas[header]+=line[:-1]
	nameHandle.close()
	return(fastas) 
