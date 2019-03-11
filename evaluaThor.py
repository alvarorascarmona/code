#! /usr/bin/env python3

##################################################


#   __________________________________
#   |                                |
#   |                                |
#   |           evaluaThor           |
#   |                                |
#   |                                |
#   |________________________________|
#               |      /|
#               |     / |
#               |    /  |
#               |   /   |
#               |  /    |
#               | /     |
#               |/      |
#               |      /|
#               |     / |
#               |    /  |
#               |   /   |
#               |  /    |
#               | /     |
#               |/      |
#               |_______|


##################################################

# Technical data:
## Author name: Álvaro Ras Carmona
## Project name: Standarization of biological measurements
## Software (code) name: evaluaThor
## Version: 1.1.0 
## Python version: Python 3.7.2
## Sysitem version (operation system used to develop evaluaThor): macOS High Sierra version 10.13.6 (17G65) 
## To run the code (2 options): ./evaluaThor.py OR python3 evaluaThor.py

############################################

'''
I tried to make the application the more interactive as posible.
That means that all along the code there is going to be a lot of prints and inputs.
By this way, we are asking all the time to the user what he/she wants to do.
'''

#################################################
#The first thing we are doing is to import the different packages
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.SubsMat.MatrixInfo import *
from Bio import SubsMat
import pandas as pd
from rpy2 import robjects
# Information obtained from: https://www.pybonacci.org/2015/06/18/trabajando-con-python-y-r/ :
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import matplotlib.pyplot as plt
from rpy2.robjects.packages import importr
##############
import os
#################################################


#################################################
###FUNCTIONS###
def ras_score (seq1,seq2,matrix_of_interest,open_penalty,extend_penalty):
    '''
    This function returns the Ras score (score invented by the autor of the script).
    It reurns a number: type = int
    What are the inputs:
        seq1 and seq2: The sequences to aling (from which score will calculate the ras score).They have to be strings.
        matrix of interest: What matrix you are going to use in the alignment. Have to be a string.
        open penalty and extend penalty that will be used in the alignment. Have to be of int type.
    Output: An int type (the score)
    '''
    alignment_1 = pairwise2.align.globalds(seq1, seq1, matrix_of_interest, open_penalty, extend_penalty)
    x1 = alignment_1[0]
    x1 = x1[2]
    alignment_2 = pairwise2.align.globalds(seq2, seq2, matrix_of_interest, open_penalty, extend_penalty)
    x2 = alignment_2[0]
    x2 = x2[2]
    alignment_problem = pairwise2.align.globalds(seq1, seq2, matrix_of_interest, open_penalty, extend_penalty)
    xf = alignment_problem [0]
    xf = xf[2]
    result = ((((xf*x2)/(x1*x2))+((xf*x1)/(x1*x2)))/2)
    result = float(result)
    return result

#Function mainly use in the "comparison_function" function
def matrix_format (number_of_matrix):
    '''
    This function returns the score matrixes you whant to compare in the apropiate format.
    It ask the user what score matrixes are going to aplicate and transform into a list.
    The only argument that needs is how many score matrixes are you going to compare (int type).
    int -> list
    '''
    matrix_comparison = ""
    for i in range (number_of_matrix):
        print ("\nAdd the matrix name (in lower case):")
        matrix_comparison += input () + ","
    matrix_comparison = matrix_comparison[0:-1]
    list_matrix = matrix_comparison.split(",")
    return list_matrix


def comparison_function (seqA,seqB,open_penalty2,extend_penalty2,number_of_matrix):
    '''
    This function is used to compare diferent score matrixes.
    It returns a data frame.
    What are the inputs:
        seq1 and seq2: The sequences to aling (from which score will calculate the ras score).They have to be strings.
        open penalty and extend penalty: Have to be of int type.
        number_of_matrix: It referes to how many score matrixes you are going to use. It has to ve a number (type int).
    Output: A dataframe with two columns (The name of the score matrix that you are comparing and the Ras score for each alignment matrix)
    '''
    matrix_use = matrix_format(number_of_matrix)
    result = ""
    for i in range(len(matrix_use)):
        matrix = matrix_use[i]
        M = mappings[matrix]
        result += str (ras_score(seqA,seqB,M,open_penalty2,extend_penalty2)) + ","
    result = result.split(",")
    result = result[0:-1]
    #Creation of the data frame
    d = {"Matrix_used" : matrix_use, "Scores" : result}
    df = pd.DataFrame(data=d)
    return df


#################################################

#################################################
###THE MAIN CODE###


print ("\nWelcome to evaluaThor software")
print ("IMPORTANT: You will need 'R' to use this software")
print("It is preferable to have 'git' installed\n")
stop_option = ""

#Here we stablish all the matrix the user can apply
mappings = {'benner22' : benner22,
'benner6' : benner6,
'benner74': benner74,
'blosum100' : blosum100,
'blosum30' : blosum30,
'blosum35' : blosum35, 
'blosum40' : blosum40, 
'blosum45' : blosum45,
'blosum50' : blosum50, 
'blosum55' : blosum55,
'blosum60' : blosum60, 
'blosum62' : blosum62,
'blosum65' : blosum65,
'blosum70' : blosum70,
'blosum75' : blosum75, 
'blosum80' : blosum80, 
'blosum85' : blosum85, 
'blosum90' : blosum90,
'blosum95' : blosum95,
'pam120' : pam120,
'pam180' : pam180,
'pam250' : pam250,
'pam30' : pam30,
'pam300': pam300,
'pam60' : pam60,
'pam90' : pam90}

#This while loop was make to avoid that the user exit from the program when he/she finish a job
while stop_option != "exit":
    #Here we are establish the main menu
    print("\nWhat do you whant to do?:")
    print ("Type 'R' to instal R software (only foy linux users)")
    print("Type 'suplement' to obtain all the suplementary material: power point, video and example .txt fasta files (only foy linux/macOS users)")
    print("Type 'use' to show how to use the software in other computer (only foy linux/macOS users)")
    print("Type 'matrix' to see all the posible matrix you can use")
    print("Type 'score' to obtain the Ras score for only 1 score matrix")
    print("Type 'compare' to compare the differents Ras scores")
    print("Type 'exit' to quit the program\n")
    main_answer = input ()

    #For the users that do not have R, we install it.
    #This option can be only used by linux users (due to the comand line for the istalation only works in linux)
    if main_answer == "R":
        os.system("sudo apt -get install r-base")
        os.system("sudo -i R")

    #We can divide this "suplement" can be divide in two parts 
    if main_answer == "suplement":
        #The firt part ask the user if he/she has git istalled
        #If not, the software install it: due to the comand line for the istalation only works in linux.
        print("\nIt is necessary to have git installed")
        print("Do you whant to install git in your system?")
        install_aswer = input("Type 'yes' (only for linux users) or 'no': ")
        if install_aswer == "yes":
            os.system("sudo apt-get install git")
            print("Git has been installed correctly")
        #The second part make a "git clone" of the repository where is all the suplementary material
        os.system("git clone https://github.com/alvarorascarmona/suplementary_material.git")
        print("\nThe suplementary material has downloaded correctly")
        print("Now you have a new directory called 'suplementary_material' with all the files")
        #Also shows all the directories that are inside the working one (you should see a new one called "suplementary_material")
        print("Your directories:")
        os.system("ls -l")
        print("\nWhat have you download:")
        os.system("ls -l suplementary_material")

    #With the "use" option the sofware shows the user the bash comand lines that has to type in the other computer where he/she whants to install it.
    if main_answer == "use":
        print("\nFor this, in the other computer it is necesary to have 'git' installed")
        print("To make the git installation apply this code (only for linux users): 'sudo apt-get install git'")
        print("To use this software in other compute follow the following code:")
        print("'git clone https://github.com/alvarorascarmona/code.git'")
        print("'cd code'")
        print("'./evaluaThor.py'")

    #With this option we can see all substitution matrix we can use 
    if main_answer == "matrix":
        posible_matrix = dir(MatrixInfo)
        #Option -> Only is saved the different score matrixes the user can use (from "dir(MatrixInfo)")
        options = posible_matrix[9:12] + posible_matrix[12:28] + posible_matrix[39:-3]
        print("\nWhat do you whant to do?")
        print("Type 'A' to see all the posible matrix")
        print("Type 'B' to see what benner, blosum or PAM matrix there are")
        print("Type 'C' to see the specific values of a score matrix")
        user_anwer = input ()
        if user_anwer == "A":
            print(options)
        if user_anwer == "B":
            print("\nExactly, in what matrix do you whant to focus?")
            user_anwer2 = input()
            if user_anwer2 == "benner":
                    print(posible_matrix [9:12])
            if user_anwer2 == "blosum":
                    print(posible_matrix [12:28])
            if user_anwer2 == "pam":
                    print(posible_matrix [39:-3])
        if user_anwer == "C":
            user_anwer3 = input("\nIndicate which one (in lower case):")
            print(SubsMat.SeqMat(mappings[user_anwer3]))

    #With this option it calculate the Ras score for one alignment
    if main_answer == "score":
        print("\nDo you have the first sequence in a .txt fasta file?")
        file_aswer = input("Type 'yes' or 'no': ")
        if file_aswer == "yes":
            #It allows to read the sequence from a fasta file (apply to al the sequence inputs)
            seq1 = ""
            file_input = input("\nIndicate the name of the file. If it is not in the same directory, indicate the path: ")
            with open(file_input) as f:
                for line in f:
                    if not line.startswith(">"):
                        seq1 += line.strip()
        if file_aswer == "no":
            #Also it allows to type by hand the sequence you whant to use (apply to al the sequence inputs)
            seq1 = input("\nAdd the first sequence (in capital letters)\n")
        print("\nDo you have the second sequence in a .txt fasta file?")
        file_aswer = input("Type 'yes' or 'no': ")
        if file_aswer == "yes":
            seq2 = ""
            file_input = input("\nIndicate the name of the file. If it is not in the same directory, indicate the path: ")
            with open(file_input) as f:
                for line in f:
                    if not line.startswith(">"):
                        seq2 += line.strip()
        if file_aswer == "no":
            seq2 = input("\nAdd the second sequence (in capital letters)\n")
        z = input("\nWhat matrix are you going to use? (in lower case): ")
        #Obtention the information from the "mapping" dictionary
        matrix_of_interest = mappings[z]
        open_penalty = int(input("\nChose the open penalty (value lower than 0): "))
        extend_penalty = int(input("\nChose the extend penalty (value lower than 0): "))
        print("\nYour Ras score is of:")
        #Aplication of the "ras_score" function
        print(ras_score(seq1,seq2,matrix_of_interest,open_penalty,extend_penalty))
        print("\nDo you whant to see the alignment/alignments?")
        answer = input("Type 'yes' or 'no': ")
        if answer == "yes":
            #Here we can see ALL the alignments performed
            alignment_problem = pairwise2.align.globalds(seq1, seq2, matrix_of_interest, open_penalty, extend_penalty)
            for i in range(len(alignment_problem)):
                print(pairwise2.format_alignment(*alignment_problem[i]))

    #With this option we obtain the different scores for the different score matrixes
    if main_answer == "compare":
        print("\nDo you have the first sequence in a .txt fasta file?")
        file_aswer = input("Type 'yes' or 'no': ")
        if file_aswer == "yes":
            seqA = ""
            file_input = input("\nIndicate the name of the file. If it is not in the same directory, indicate the path: ")
            with open(file_input) as f:
                for line in f:
                    if not line.startswith(">"):
                        seqA += line.strip()
        if file_aswer == "no":
            seqA = input("\nAdd the first sequence (in capital letters)\n")
        print("\nDo you have the second sequence in a .txt fasta file?")
        file_aswer = input("Type 'yes' or 'no': ")
        if file_aswer == "yes":
            seqB = ""
            file_input = input("\nIndicate the name of the file. If it is not in the same directory, indicate the path: ")
            with open(file_input) as f:
                for line in f:
                    if not line.startswith(">"):
                        seqB += line.strip()
        if file_aswer == "no":
            seqB = input("\nAdd the second sequence (in capital letters)\n")
        open_penalty2 = int(input("\nChose the open penalty (value lower than 0): "))
        extend_penalty2 = int(input("\nChose the extend penalty (value lower than 0): "))
        number_of_matrix = int(input("\nHow many alignment matrix do you whant to compare?: "))
        x = comparison_function(seqA,seqB,open_penalty2,extend_penalty2,number_of_matrix)
        print("\nHere is your result:")
        print(x)
        print("\nDo you whant to see the statistics of the scores?")
        see_statistics = input("Type 'yes' or 'no': ")
        #Here we create a .txt with the output of the "comparison_function"
        #This .txt will be use if the R code is run
        with open("tabla_scores.txt", "w") as y:
            y.write(str(x))
        
        #Here, with R code, we can obtain the statistics of the dataframe 
        if see_statistics == "yes":
            codigo_summary = """
            wd <- getwd()
            setwd(wd)
            x <- read.table ("tabla_scores.txt", header = T)
            summary(x)
            """
            print(ro.r(codigo_summary))

        #Here, with R code, we plot that data frame (using a bar plot)
        print("\nDo you whant to see the plot of the scores?")
        see_plot = input("Type 'yes' or 'no': ")
        if  see_plot == "yes":
            codigo_plot = """
            wd <- getwd()
            setwd(wd)
            x <- read.table ("tabla_scores.txt", header = T)
            pdf("Scores_histogram.pdf")
            barplot(x$Scores, names = x$Matrix_used, col = "light blue", xlab = "Matrix used", ylab = "Score", main = "Score comparison")
            dev.off()
            """
            print(ro.r(codigo_plot))
            print("\nYour plot has been correctly generated.")
            print("Please, check the your working directory.")
            print("You should have a PDF called 'Scores_histogram.pdf'.")

        #The .txt generated can be saved in a .txt file
        print("\nDo you whant to save the results (score matrixes used vs Ras scores) on a .txt file?")
        save_option = input("Type 'yes' or 'no': ")
        if save_option == "no":
            #With this code lines we delete the .txt generated before
            path = os.getcwd() 
            y = path + "/" + "tabla_scores.txt"
            os.remove(y)
        else: 
            print("\nThe file has been corretly saved with the 'tabla_scores.txt' name")
            print("Please, check your working directory")

    #Here we let the user to get out of the program (get out of the while loop)
    if main_answer == "exit":
        print("\nAre you sure you want to exit?")
        print("Type 'yes' or 'no'\n ")
        second_oportnity = input()
        if second_oportnity == "yes":
            stop_option = "exit"


print("\nGoodbye")

      
#################################################