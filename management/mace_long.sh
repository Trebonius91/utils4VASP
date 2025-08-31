#!/bin/sh
#
#    mace_long.log: This script manages long MACE MD simulations with
#                   several restarts on a cluster with walltime limit.
#
#    Author: Maximilian Bechtel  2025 (maxi.bechtel@fau.de)
#
#    Part of utils4VASP
#
#    Adaption of md_long.log for MACE MD simulations
#    by Julien Steffen           2023 (julien.steffen@fau.de)
#    by Andreas Moelkner         2023 (andreas.moelkner@fau.de)
#
#
#   NB: the python script modify_poscar.py -cart2frac is used
#   to generate the new POSCAR file in fractional coordinates
#   simply 'mv CONTCAR POSCAR' will slow down the calculation
#   since the cartesian coordinates will grow bigger and bigger
#   however: this will result in new initial velocities for each
#            restart, this might be changed if needed
#
#   to ensure that the script can run in the background on a cluster
#   use tmux and load a python environment before starting the script
#   e.g. by 'module load python' in the user´s .profile file

##########################################################
### GLOBAL VARIABLES
##########################################################

# Check every $sleep_time minutes the status of the calcuation
sleep_time=300
# Counter for the MD steps already performed after each abortion
nsteps_sum=0
# Check if calculation has performed all desired MD steps
finished=0

#
### Usage of the script
#
echo "This script manages long MACE MD simulations with several restarts
on a cluster with walltime limit."
echo ""
echo "- If a calculation is aborted due to the time limit, the current
  CONTCAR, XDATCAR, POSCAR, mace.log, dummy, gradient.dat and dynamic.py
  file will be renamed with indices of the respective part of the calculation."
echo "- If the calculation has been completed, the script will exit
  normally."
echo "- The desired maximum number of restarts must be limited as command
  line argument."
echo "- If the script is started on a cluster, the best would be to run it in
  a seperate tmux window. Otherwise, the script will be cancel automatically
  when logging out of the cluster terminal."
echo "- If the script shall be terminated before the maximum number of restarts
  are reached, a file 'kill' must be created by 'touch kill'."
echo "- The current status of the calculation will be written to the file
  'md_long.log'"
echo "- The python script modify_poscar.py -cart2frac is used generate the new
  POSCAR file in fractional coordinates, simply 'mv CONTCAR POSCAR' will slow
  down the calculation since the cartesian coordinates will grow bigger and bigger
  however: this will result in new initial velocities for each restart, this might
  be changed if needed"
echo "- To ensure that the script can run in the background on a cluster, use tmux
  and load a python environment before starting the script, e.g. by 'module load python'
  in the user´s .profile file"

echo ""
echo "Usage: md_long.sh [max no. of restarts]"
echo ""

#
### Check if max. no. of restartes are given
### by checking if string "$1" is empty
#
if [ -z "$1" ]; then
  echo "Max. no. of restarts was not given!"
  echo ""
  exit 1
fi
num_restarts=$1

#
### remove file md_long.sh from a previous run if it exits
#
if test -f "md_long.log"; then
  rm md_long.log
  echo "Remove old file md_long.log ..."
fi
#
### remove kill file if user forgot to delete an old one
#
if test -f "kill"; then
  rm kill
  echo "Remove old kill file ..."
fi

touch md_long.log
echo "This is the logfile of mace_long.sh." > md_long.log
echo "The calculation will be restarted $num_restarts times until the calculation will be finished." >> md_long.log
echo "Generate a file 'kill' if the script shall be aborted." >> md_long.log
echo " --- " >> md_long.log

#
### Get the desired total number of MD steps
### nsteps in dynamic.py
#
string=$(sed -nr "/nsteps /p" dynamic.py)
string_array=( $string )
nsteps_total=$( echo ${string_array[2]} )
echo "Found total number of MD steps: $nsteps_total"
echo "Total number of MD steps: $nsteps_total" >> md_long.log

##########################################################
### THE MAIN LOOP
##########################################################
echo " "
echo "Entering main loop ..."
echo ""
echo " " >> md_long.log
echo "Entering main loop ..." >> md_long.log
echo "" >> md_long.log

#
### OUTER LOOP: restart calculation if needed
### but only $num_restarts times
#
for ((i=1;i<=$num_restarts;i++))
do
  #
  ### Submit next calculation run and get the current jobnumber
  #
  sbatch slurm_script > submit.txt
  jobnum=$(awk '{print $4}' submit.txt)
  echo "Job number $jobnum started." >> md_long.log
  rm submit.txt

  #
  ### INNER LOOP: one calculation run
  ### break inner loop if walltime limit reached
  #
  minute=0
  for (( ; ; ))
  do
    #
    ### Check if user created the kill file
    #
    if test -f "kill"; then
      echo "The 'kill' file has been found."
      echo "The 'kill' file has been found." >> md_long.log
      echo "Job number $jobnum will be canceled ..."
      echo "Job number $jobnum will be canceled ..." >> md_long.log
      scancel $jobnum
      echo "mace_long.log will be terminated ..."
      echo "mace_long.log will be terminated ..." >> md_long.log

      exit 0
    fi

    #
    ### Check if calculation already finished
    #

    if grep -Fq "The calculation has been finished successfully!" mace.log
    then
      finished=1
    fi

    #
    ### Check if calculation still exists in queue
    ### if not break inner loop
    #

    # Try to get the jobnumber via squeue
    jobnum_found=$(squeue | grep $jobnum | awk '{print $1}')
    # Get status of job
    job_status=$(squeue | grep $jobnum | awk '{print $5}')
    # Get the runtime of the job
    job_runtime=$(squeue | grep $jobnum | awk '{print $6}')

    if [ "$jobnum_found" == "$jobnum" ]  # check if the correct jobnumber is found
    then

      if [ "$job_status" == "R" ]  # check if the job is running
      then
        echo "Run $i of $num_restarts, minute $minute: $jobnum has status $job_status and Runtime $job_runtime" >> md_long.log
      elif [ "$job_status" == "PD" ]  # check if the job is still pending 
      then
        echo "Run $i of $num_restarts, minute $minute: $jobnum has status $job_status" >> md_long.log
      fi

    else
      #
      ### Wait for another $sleep_time minutes, maybe the squeue command is currently unavailable
      #
      echo "Run $i of $num_restarts, minute $minute: Job $jobnum not found." >> md_long.log
      echo "Waiting for another 5 minutes and check again ..." >> md_long.log

      sleep $sleep_time
      minute=$(($minute + 5))

      if [ "$jobnum_found" == "$jobnum" ]
      then
        echo "Run $i of $num_restarts, minute $minute: Job $jobnum was found again." >> md_long.log
      else
        echo "Run $i of $num_restarts, minute $minute: The calculation will be restarted." >> md_long.log
        break  # break the inner loop to restart the calculation
      fi

    fi

    #
    ### Wait for $sleep_time minutes and increase the minute counter
    #
    sleep $sleep_time
    minute=$(($minute + 5))
  done

  #
  ### Check if calculation was finished within the last run
  #

  if [ $finished -eq 1 ]
  then
    echo "HURRAY! Your MD simulation has been finished!" >> md_long.log
    break  # break the outer loop since the calculation has been finished already
  fi

  #
  ### Get the current number of MD steps when calculation aborted
  ### due to the time limit
  #
  string=$(tail -1 mace.log)
  string_array=( $string )
  nsteps_current=$(echo ${string_array[0]})
  nsteps_sum=$(($nsteps_sum + $nsteps_current))

  echo "-------Run $i of $num_restarts finished! ------------------" >> md_long.log
  echo "MD steps finished: $nsteps_sum of $nsteps_total" >> md_long.log
  echo "" >> md_long.log

  #
  ### Prepare the new input files
  #
  cp dynamic.py dynamic$i.py
  cp mace.log mace$i.log
  cp dummy dummy$i
  cp gradients.dat gradients$i.dat
  cp XDATCAR XDATCAR$i

  #
  ### Translate the current CONTCAR to the new POSCAR file
  ### NB: the python script modify_poscar.py is used
  ### to generate the new POSCAR in fractional coordinates
  ### simply moving CONTCAR to POSCAR will slow down the calculation
  ### since the cartesian coordinates will grow bigger and bigger
  #
  cp POSCAR POSCAR$i
  cp CONTCAR CONTCAR$i
  mv CONTCAR POSCAR
  modify_poscar.py -cart2frac
  mv POSCAR_mod POSCAR

  #
  ### Update the number of MD steps in dynamic.py
  #  
  echo "Updating the number of MD steps in dynamic.py ..."
  nsteps_new=$(($nsteps_total - $nsteps_sum))

  line_number=$(sed -n /nsteps/= dynamic.py | head -1)  # line number of nsteps
  # Insert new MD steps after $line_number
  sed -i "$line_number a\nsteps = $nsteps_new" dynamic.py
  # Delete old MD steps on $line_number
  sed -i "$line_number d" dynamic.py

done

echo ""
echo "The script mace_long.sh has finished! Goodbye!" >> md_long.log
