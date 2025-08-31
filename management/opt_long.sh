#!/bin/sh
#
#    opt_long: Management of long VASP geometry 
#      optimizations with several restarts on a server
#      with wall time restrictions
#      During each restart, the input for the next calculation
#      (INCAR, POSCAR, ML_AB) is built from the output of the 
#      last calculation. Further, several errors or unexpected
#      behaviors can be managed
#    Part of VASP4CLINT
#      Julien Steffen, 2023 (julien.steffen@fau.de)
#

echo "
 This script manages long VASP geometry optimization 
  calculations with several restarts. 
 If a calculation is stopped due to time limit,
  the current CONTCAR file will be copied to a new POSCAR
  file, and OUTCAR and XDATCAR will be renamed with indices 
  of the respective parts of the calculation.
 If the calculation has been completed, the script will exit normally.
 The desired maximum number of restarts must be given by the
  user as a command line argument:
 Generate file 'kill' to terminate the script ahead of schedule.
 The current status of the script will be written to 'opt_long.log'

 Usage: opt_long.sh [No. of restarts]
"
if [ -z "$1" ]; then
   echo "No command line argument submitted!"
   exit 1
fi
restarts=$1


if test -f "opt_long.log"; then
   rm opt_long.log
fi
	
touch opt_long.log

#
#    Save initial INCAR file (since the actual one will be changed)
#    Also save initial POSCAR file
#
cp INCAR INCAR_initial
cp POSCAR POSCAR_initial

echo "This is the logfile of 'md_long.sh'" > opt_long.log
echo "$restarts runs will be started, each until calculation aborted or finished" >> opt_long.log
echo "Generate a file 'kill' if the script shall be aborted." >> opt_long.log
echo " --- " >> opt_long.log

#Define sleep time and minute counter
sleept=120 #Multiples of 60!!!!
min=$(($sleept / 60 ))


istep_sum=0
finished=0

#Define the different possible error messages
re1="No space left on device"
re2="srun: Job step aborted:"
re3="NaN"
re4="Firewall refused connection."
re5="se ML_MB !!!"

rm vasp.err
rm vasp.out
for ((i=1; i<=$restarts;i++))
do

   sbatch slurm_script > submit.txt
   jobnum=$(awk '{print $4}' submit.txt)
   echo "Job number $jobnum started." >> opt_long.log

   minute_act=0
   for (( ; ; ))
   do 
 #
 #     Stay in loop for current calculation until time limit is reached or calculation
 #     is finished!
 #  
      echo "Run $i of $restarts, Minute $minute_act " >> opt_long.log
      sleep 60
#    Calculation aborted 
      if grep -Fq "DUE TO TIME LIMIT" vasp.err
      then
	 rm vasp.err
         break
      fi
#    Calculation finished
      if grep -Fq "General timing and accounting " OUTCAR
      then
         rm vasp.err
	 finished=1
         break
      fi
#    Check for several different VASP errors, restart current calculation if needed 

      message=$(cat vasp.err)
      message2=$(tail -1 vasp.out)
   # If the error occurred, restart the calculation
      if [[ $message =~ $re1 ]]
      then
	 echo "Job No. $jobnum will be canceled."  >> opt_long.log
         rm vasp.err 
         rm vasp.out
         rm vasprun.xml
    #     rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'No space left on device' error!" >> opt_long.log
	 echo "Job number $jobnum started." >> opt_long.log
      fi

      if [[ $message =~ $re2 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> opt_long.log
         rm vasp.err
         rm vasp.out
         rm vasprun.xml
   #      rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'srun: Job step aborted' error!" >> opt_long.log
         echo "Job number $jobnum started." >> opt_long.log
      fi


      if [[ $message =~ $re4 ]]
      then
	 echo "Job No. $jobnum will be canceled."  >> opt_long.log
         rm vasp.err 
         rm vasp.out
         rm vasprun.xml
   #      rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'Firewall refused connection.' error!" >> opt_long.log
	 echo "Job number $jobnum started." >> opt_long.log
      fi

      # If the NaN error occured, abort and restart the calculation
      if [[ $message2 =~ $re3 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> opt_long.log
         scancel $jobnum
         rm vasprun.xml
   #      rm OUTCAR
         rm vasp.err
         rm vasp.out
         mv XDATCAR XDATCAR_part$i
         sleep $sleept
         sbatch slurm_script > submit.txt
	 jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'NaN position' error!" >> opt_long.log
	 echo "Job number $jobnum started." >> opt_long.log
      fi

#      
#Checking in squeue for the status of the calculation with jobnum
#
      job_number=$(squeue | grep $jobnum | awk '{print $1}') #get jobnumber
      job_status=$(squeue | grep $jobnum | awk '{print $5}') #get status
      job_runtime=$(squeue | grep $jobnum | awk '{print $6}') #get runtime


      if [ $jobnum -eq $job_number ] #check if jobnum is found
      then
         if [ "$job_status" == "R" ] #check if running
         then
            echo "Job $jobnum has status $job_status and Runtime $job_runtime" >> opt_long.log
         elif [ "$job_status" == "PD" ] #check if pending
         then
            echo "Job $jobnum has status $job_status" >> opt_long.log
         else
            if [ -z "$jobnum" ]
            then
               echo "Variable $jobnum is empty! The job will be restarted..." >> opt_long.log
               cp vasp.err vasp.err$i.$minute_act
               cp vasp.out vasp.out$i.$minute_act
               rm vasprum.xml
               rm OUTCAR
               rm vasp.err
               rm vasp.out
               sleep $sleept
               sbatch slurm_script > submit.txt
               jobnum=$( awk '{print $4}' submit.txt)
               echo "Job number $jobnum started." >> opt_long.log
            else
               echo "Job $jobnum is getting cancled (or other weird status): Status $job_status" >> opt_long.log
            fi
         fi
      else
         echo "Job $jobnum not found" >> ml_long.log

         #Restart with current input
         echo "Restarted due to unknown abortion!" >> opt_long.log

         cp vasp.err vasp.err$i.$minute_act
         cp vasp.out vasp.out$i.$minute_act
         rm vasprum.xml
         rm OUTCAR
         rm vasp.err
         rm vasp.out
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$( awk '{print $4}' submit.txt)

         echo "Job number $jobnum started." >> opt_long.log
      fi



      if test -f "kill"; then
         echo "The 'kill' file has been found! "
	 echo " md_long.sh will be terminated..."
	 echo "The 'kill' file has been found! " >> opt_long.log
	 echo " md_long.sh will be terminated..." >> opt_long.log
	 scancel $jobnum
	 exit 0
      fi
      minute_act=$(($minute_act + $min))
   done	   

   if [ $finished -eq 1 ]
   then
      echo "HURRAY! Your MD simulation has been finished!" >> opt_long.log
      break
   fi	   
  
   echo "-------Run $i of $restarts finished! ------------------" >> opt_long.log
   cp CONTCAR CONTCAR$i
   cp CONTCAR POSCAR
   cp OUTCAR OUTCAR$i
   cp XDATCAR XDATCAR$i
   
   sleep $sleept
done	

echo "  " >> opt_long.log
echo " The script 'opt_long.sh' has finished! Goodbye! " >> opt_long.log

