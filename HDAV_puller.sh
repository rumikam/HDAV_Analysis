#!/bin/bash
#
#SBATCH --job-name=HDAV_puller                                                                                               #Job name
#SBATCH --mail-type=ALL   
#SBATCH --mail-user=rumika.mascarenhas@ucalgary.ca                                                                                    #Where to send mail
#SBATCH --ntasks=1                                                                                                              #Run on single CPU
#SBATCH --partition=apophis
#SBATCH --nodes=1            #Job memory request
#SBATCH --time=6:00:00                                    #wall time limit
#SBATCH --output=serial_HDAV_%j.log                                             #Standard output/error logs


# this is where you chage your directory to where your script is
# give that your script already writes to file, you dont need to change the output flow
module load bioconda/2018.11 
python HDAV_puller.py


#REALLY IMPORTANT
#To submit a slurm job, you need to go sbatch myscript.sh
#To check on a job type squeue --> this will show all jobs running on the cluster
