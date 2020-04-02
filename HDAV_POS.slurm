#!/bin/bash
#SBATCH --job-name=Pre-process
#SBATCH --mail-type=ALL                                                                                             #Job name
#SBATCH --mail-user=rumika.mascarenhas@ucalgary.ca                                                                                    #Where to send mail
#SBATCH --ntasks=1                                                                                                          #Run on single CPU
#SBATCH --partition=apophis
#SBATCH --nodes=1    
#SBATCH --time=10:00:00      #wall time limit
#SBATCH --output=pre-process_%j.log

python HDAV_position.py
