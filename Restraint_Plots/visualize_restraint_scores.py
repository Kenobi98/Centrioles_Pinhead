#############$Plot the score values for Restraints in stat file$#############
#@Author: Kartik Majila
#......Use the shell script "read_stat_out.sh" to run this python script recursively

####@One Above All@####
import numpy
import matplotlib.pyplot as plt
import sys

file = sys.argv[1]                 #Score file created using process_output.py
name = sys.argv[2]                 #Name of the plot image

print("Preparing to plot!")

#Open the file in reading mode
with open (file,"r") as f:
    score_file = f.readlines()

score = []
score_final = []

#Read the score values from the score file created
[score.append(score_file[i].split(" ")[1].strip()) for i in range(len(score_file)) ]

#Remove the extra '' elements in the list
if name == "GaussianEMRestraint_None": 
    [score_final.append(float(score[i])) for i in range(0,len(score), 3)]
else:
    [score_final.append(float(score[i])) for i in range(0,len(score), 2)]

frames = range(len(score_final))

#Plot the score vs no. of frames
plt.scatter(frames, score_final)
plt.title(name)
plt.xlabel("No. of Frames")
plt.ylabel("Score ")
plt.savefig(name+".png")

print("Plot ready for_", name)