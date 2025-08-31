#!/bin/sh
#
#    ml_long: Management of long VASP Machine learning
#      calculations with several restarts on a server 
#      with wall time restrictions. 
#      During each restart, the input for the next calculation
#      (INCAR, POSCAR, ML_AB) is built from the output of the 
#      last calculation. Further, several errors or unexpected
#      behaviors can be managed
#    Part of VASP4CLINT
#      Julien Steffen, 2023 (julien.steffen@fau.de)
#      Andreas Mölkner, 2023 (andreas.moelkner@fau.de)
#

echo "
 This script manages long VASP Machine Learning calculations 
   with several restarts.
 If a calculation is stopped due to time limit,
   the current CONTCAR file will be copied to a new POSCAR
   file, and OUTCAR, XDATCAR and ML_LOGFILE will be renamed with indices 
   of the respective parts of the calculation.
 If the calculation has been completed, the script will exit normally.
   The parameter ML_CTIFOR will be set based on the output of 
   the previous run in order to get a smooth learning curve.
 Machine learnings under influence of applied forces (steered MD) can
   be handled as well. If steered MD learnings shall be done, the keywords
   SPRING_K, SPRING_R0 and SPRING_V0 must be added, together with a 
   suitable ICONST file (see VASP wiki=.  
 The start and end temperatures will be set based on the number 
   of already calculated steps and the global TEBEG and TEEND
 The desired maximum number of restarts must be given by the
   user as a command line argument:
 INFO: CONTCAR files are now translated to POSCAR files with 
   the contcar2poscar program. If problems appear, this might be 
   changed! 
 Generate file 'kill' to terminate the script ahead of schedule.
   The current status of the script will be written to 'ml_long.log'
 
 Usage: ml_long.sh [No. of restarts] [Steered Option]
 
 Possible options (referring to steered MD) are: 
   '-no_steered'  No bias forces will be applied
   '-volume_steered'  A bias force applied on the volume 
            of the simulation cell will be considered. All three 
	    coordinate axis must be given in the SPRING... commands
	    Further, angles should be constrained in ICONST file
   '-atom_steered' A bias force will be applied on one or several
            coordinates of one of several atoms 
 "

# Input of Restarts

if [ -z "$1" ]; then
   echo "No restart command line argument submitted!"
   exit 1
fi
restarts=$1

# Option for steered MD

if [ -z "$2" ]; then
   echo "No steered MD command line argument submitted!"
   exit 1
fi
steered_opt=$2

if [ $steered_opt == "-no_steered" ]; then
   echo "No steered MD will be applied."
elif [ $steered_opt == "-volume_steered" ]; then  
   echo "The volume of the cell will be changed during learning."
elif [ $steered_opt == "-atom_steered" ]; then
   echo "One or several atoms will be steered during learning."
else 
   echo "No valid option for steered MD was given!"
   exit 1
fi	


#Define sleep time and minute counter 
sleept=120 #Multiples of 60!!!!
min=$(($sleept / 60 ))

####################################################################################################################################
#################### FUNCTIONS #####################################################################################################


#define the continue function -> Adapting INCAR and other files to continue a calculation after abortion/error:
function con_calc { 
   string=`sed -n '/T=/h; ${x;p;}' vasp.out`
   string2=( $string )
   step_act=`echo ${string2[0]}`
   step_sum=$(($step_sum + $step_act))

   tebeg_act=$(($tebeg + ($teend - $tebeg) * $step_sum/$step_all))

   nsw_act=$(($step_all - $step_sum))


   cp CONTCAR CONTCAR$i.$minute_act
   cp CONTCAR POSCAR
   contcar2poscar
   cp OUTCAR OUTCAR$i.$minute_act
   cp XDATCAR XDATCAR$i.$minute_act
   cp ML_LOGFILE ML_LOGFILE$i.$minute_act
   cp ML_ABN ML_ABN$i.$minute_act
   cp ML_ABN ML_AB


   #  Remove old ML_CTIFOR and ML_ISTART parameters, if present
   sed -i '/ML_CTIFOR/d' ./INCAR
   sed -i '/ML_ISTART/d' ./INCAR

   # Remove old TEBEG parameter 

   sed -i '/TEBEG/d' ./INCAR

   # Remove old NSW (step number) parameter
   sed -i '/NSW/d' ./INCAR

   # update the ML_CTIFOR parameter for the next trajectory part 
   ctifor_old=$ctifor
   ctifor=` grep "THRUPD" ML_LOGFILE | tail -1 | awk '{ printf( $4) }'`
   # If the variable is empty, use the previous value 
   if [ -z "$ctifor" ]; then
      ctifor=$ctifor_old
   fi
   echo "ML_CTIFOR = $ctifor" >> INCAR
   #echo "ML_ISTART = 1" >> INCAR

   # Add new starting temperature
   echo "TEBEG = $tebeg_act" >> INCAR

   # Add new MD step number
   echo "NSW = $nsw_act" >> INCAR

}

function update_ml_mconf {

   if test -f "ML_ABN" ; then
      nconf=$(head -n +5 ML_ABN | tail -n -1) #reads number of conformers out of ML_ABN
      inconf=$(grep 'ML_MCONF' INCAR | awk '{print $3}') #gets ML_MCONF
      diff=$(($inconf - $nconf)) #calculates difference
      outconf=$(($inconf + 1000))

      if [ "$diff" -lt 1000 ]
      then
         echo "Increase ML_MCONF from $inconf to $outconf"
         sed -i '/ML_MCONF/d' ./INCAR
         echo "ML_MCONF = $outconf" >> INCAR
     fi
   fi

}





####################################################################################################################################
############################## THE CODE BEGINNS ####################################################################################
if test -f "ml_long.log"; then
   rm ml_long.log
fi
	
touch ml_long.log

#
#    Save initial INCAR and POSCAR file (since the actual one will be changed)
#
cp INCAR INCAR_initial
cp POSCAR POSCAR_initial #in case of early restarts

#
#    Read start and end temperature from INCAR file 
#

string=`sed -nr '/TEBEG/p' INCAR`
string2=( $string )
tebeg=`echo ${string2[2]}`

string=`sed -nr '/TEEND/p' INCAR`
string2=( $string )
teend=`echo ${string2[2]}`

string=`sed -nr '/NSW/p' INCAR`
string2=( $string )
step_all=`echo ${string2[2]}`

echo "This is the logfile of 'ml_long.sh'" > ml_long.log
echo "$restarts runs will be started, each until calculation aborted or finished" >> ml_long.log
echo "Generate a file 'kill' if the script shall be aborted." >> ml_long.log
echo " --- " >> ml_long.log
echo " * Total number of MD timesteps: $step_all " >> ml_long.log
echo " * Start temperature: $tebeg K " >> ml_long.log
echo " * End temperiature: $teend K " >> ml_long.log 
echo "  " >> ml_long.log


#  Remove old ML_CTIFOR and ML_ISTART parameters, if present
#sed -i '/ML_CTIFOR/d' ./INCAR
#sed -i '/ML_ISTART/d' ./INCAR  USE ML_MODE = train

ctifor_old=0.0020000
ctifor=0.00200000
istep_sum=0
finished=0

re1="No space left on device"
re2="srun: Job step aborted:"
re3="NaN"
re4="Firewall refused connection."
re5="se ML_MB !!!"
re6="application called MPI_Abort"
re7="severe (24): end-of-file during read, unit 15"


rm vasp.err

#################################################################################################################################################################################
########### MAIN LOOP (NUMBER OF RESTARTS) ############


for ((i=1; i<=$restarts;i++))
do

#   if [ $i -eq 1 ]; then
#      echo "ML_ISTART = 0" >> INCAR
#   fi	   

   sbatch slurm_script > submit.txt
   jobnum=$(awk '{print $4}' submit.txt)
   echo "Job number $jobnum started." >> ml_long.log



   minute_act=0
   for (( ; ; )) #Inner Loop to check every Minute the state of the calculation
   do 
 #
 #     Stay in loop for current calculation until time limit is reached or calculation
 #     is finished!
 #  
      echo "Run $i of $restarts, Minute $minute_act " >> ml_long.log
      sleep $sleept
#    Calculation aborted 
      if grep -Fq "DUE TO TIME LIMIT" vasp.err
      then
	 rm vasp.err
         break
      fi
#    Calculation finished
      if grep -Fq "TIMING" ML_LOGFILE
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
	 echo "Job No. $jobnum will be canceled."  >> ml_long.log
         rm vasp.err 
         rm vasp.out
         rm vasprun.xml
         rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'No space left on device' error!" >> ml_long.log
	 echo "Job number $jobnum started." >> ml_long.log
      fi

      if [[ $message =~ $re2 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> ml_long.log


         rm vasp.err
         rm vasp.out
         rm vasprun.xml
         rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'srun: Job step aborted' error!" >> ml_long.log
         echo "Job number $jobnum started." >> ml_long.log
      fi


      if [[ $message =~ $re4 ]]
      then
	 echo "Job No. $jobnum will be canceled."  >> ml_long.log
         rm vasp.err 
         rm vasp.out
         rm vasprun.xml
         rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'Firewall refused connection.' error!" >> ml_long.log
	 echo "Job number $jobnum started." >> ml_long.log
      fi

      # If the NaN error occured, abort and restart the calculation
      if [[ $message2 =~ $re3 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> ml_long.log
         scancel $jobnum
         rm vasprun.xml
         rm OUTCAR
         rm vasp.err
         rm vasp.out
         mv XDATCAR XDATCAR_part$i
         sleep $sleept
         sbatch slurm_script > submit.txt
	 jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'NaN position' error!" >> ml_long.log
	 echo "Job number $jobnum started." >> ml_long.log
      fi

      if [[ $message =~ $re6 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> ml_long.log

         con_calc
	 update_ml_mconf

         rm vasp.err
         rm vasp.out
         rm vasprun.xml
         rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
	 echo "Restarted due to 'MPI_Abort' error! (CONTINUED)" >> ml_long.log
         echo "Job number $jobnum started." >> ml_long.log
      fi


      # If the ML_MB error occurs, restart the calculation
      if [[ $message2 =~ $re5 ]]
      then
	 string=`sed -nr '/ML_MB/p' INCAR`
         string2=( $string )
         ml_mb=`echo ${string2[2]}`
         ml_mb=$(($ml_mb + 800))
	
	 sed -i '/ML_MB/d' ./INCAR 
         echo "ML_MB = $ml_mb" >> INCAR

         rm vasprun.xml
         rm OUTCAR
         rm vasp.err
         rm vasp.out
         mv XDATCAR XDATCAR_part$i
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'increase ML_MB' error!" >> ml_long.log
         echo "Job number $jobnum started." >> ml_long.log
      fi
      
      if [[ $message2 =~ $re7 ]]
      then
         echo "Job No. $jobnum will be canceled."  >> ml_long.log
         
	 cp CONTCAR_scratch POSCAR

         rm vasprun.xml
         rm OUTCAR
         rm vasp.err
         rm vasp.out
         mv XDATCAR XDATCAR_part$i
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
         echo "Restarted due to 'end-of-file' error!" >> ml_long.log
         echo "Job number $jobnum started." >> ml_long.log
      fi
#CHEKCING IF ATOMS ARE IN THE GASPHASE

      check_vac.py
      

      if test -f "gas_kill"; then
         rm gas_kill
         echo "Atoms in the Gasphase detected " >> ml_long.log
         echo "The cacluation will be restartet with CONTCAR_scratch" >> ml_long.log

         scancel $jobnum

	 cp CONTCAR dummy 
         cp CONTCAR_scratch CONTCAR
         
         sleep 60

	 con_calc
	 update_ml_mconf

	 mv dummy CONTCAR$i.$minute_act
	 cp vasp.out vasp.out$i.$minute_act.gas
         rm vasp.err
         rm vasp.out
         rm vasprun.xml
         rm OUTCAR
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$(awk '{print $4}' submit.txt)
	 echo "Restarted due to Atom in Gasphase (CONTINUED)" >> ml_long.log
         echo "MD steps finished:  $step_sum   of $step_all" >> ml_long.log
         echo "Job number $jobnum started." >> ml_long.log

      fi






#Checking in squeue for the status of the calculation with jobnum

      job_number=$(squeue | grep $jobnum | awk '{print $1}') #get jobnumber
      job_status=$(squeue | grep $jobnum | awk '{print $5}') #get status
      job_runtime=$(squeue | grep $jobnum | awk '{print $6}') #get runtime
      
      
      if [ $jobnum -eq $job_number ] #check if jobnum is found
      then
         if [ "$job_status" == "R" ] #check if running
         then
            echo "Job $jobnum has status $job_status and Runtime $job_runtime" >> ml_long.log
         elif [ "$job_status" == "PD" ] #check if pending
         then
            echo "Job $jobnum has status $job_status" >> ml_long.log
         else
            if [ -z "$jobnum" ]
            then 
               echo "Variable $jobnum is empty! The job will be restarted..." >> ml_long.log
               cp vasp.err vasp.err$i.$minute_act
               cp vasp.out vasp.out$i.$minute_act
               rm vasprum.xml
               rm OUTCAR
               rm vasp.err
               rm vasp.out
               sleep $sleept
               sbatch slurm_script > submit.txt
               jobnum=$( awk '{print $4}' submit.txt)
               echo "Job number $jobnum started." >> ml_long.log
            else     
               echo "Job $jobnum is getting cancled (or other weird status): Status $job_status" >> ml_long.log
            fi
         fi
      else
         echo "Job $jobnum not found" >> ml_long.log

	 #Restart with current input
         echo "Restarted due to unknown abortion!" >> ml_long.log

	 cp vasp.err vasp.err$i.$minute_act
	 cp vasp.out vasp.out$i.$minute_act
         rm vasprum.xml
         rm OUTCAR
         rm vasp.err
         rm vasp.out
         sleep $sleept
         sbatch slurm_script > submit.txt
         jobnum=$( awk '{print $4}' submit.txt)

	 echo "Job number $jobnum started." >> ml_long.log
      fi

      

#Check for the KILL file and terminate ml_long and running calculation!
      if test -f "kill"; then
         echo "The 'kill' file has been found! "
	 echo " ml_long.sh will be terminated..."
	 echo "The 'kill' file has been found! " >> ml_long.log
	 echo "Job No. $jobnum will be canceled."  >> ml_long.log
	 scancel $jobnum
	 echo " ml_long.sh will be terminated..." >> ml_long.log
	 exit 0
      fi
      minute_act=$(($minute_act + $min))
   done	   #Ends inner minute Loop



################
#This is outer (Restarts) Loop again


   if [ $finished -eq 1 ]
   then
      echo "HURRAY! Your ML-FF generation has been finished!" >> ml_long.log
      break
   fi	   
#
#    Find out number of calculated time steps 
#
   string=`sed -n '/T=/h; ${x;p;}' vasp.out`
   string2=( $string )
   step_act=`echo ${string2[0]}`
   step_sum=$(($step_sum + $step_act))

   tebeg_act=$(($tebeg + ($teend - $tebeg) * $step_sum/$step_all))

   nsw_act=$(($step_all - $step_sum))
#
#    Write new TEBEG based on fraction of calculated time steps
#

  
   echo "-------Run $i of $restarts finished! ------------------" >> ml_long.log
   echo "  MD steps finished:  $step_sum   of $step_all " >> ml_long.log
   cp CONTCAR CONTCAR$i
   cp CONTCAR POSCAR
   contcar2poscar
   cp OUTCAR OUTCAR$i
   cp XDATCAR XDATCAR$i
   cp ML_LOGFILE ML_LOGFILE$i
   cp ML_ABN ML_ABN$i
   cp ML_ABN ML_AB

   update_ml_mconf
#
#     For simulations with variable cell volume (SPRING_K is present)
#     The complete cell topology is considered in calculating the new 
#     initial axis lengths 
#     Angles must be hold constant in ICONST file
#
   if [ $steered_opt == "-volume_steered" ]; then
      length_x1=`sed -n '3p' CONTCAR | awk '{ printf( $1) }'`
      length_x2=`sed -n '3p' CONTCAR | awk '{ printf( $2) }'`
      length_x3=`sed -n '3p' CONTCAR | awk '{ printf( $3) }'`
      length_y1=`sed -n '4p' CONTCAR | awk '{ printf( $1) }'`
      length_y2=`sed -n '4p' CONTCAR | awk '{ printf( $2) }'`
      length_y3=`sed -n '4p' CONTCAR | awk '{ printf( $3) }'`
      length_z1=`sed -n '5p' CONTCAR | awk '{ printf( $1) }'`
      length_z2=`sed -n '5p' CONTCAR | awk '{ printf( $2) }'`
      length_z3=`sed -n '5p' CONTCAR | awk '{ printf( $3) }'`
      length_x=$(awk -v x1=$length_x1 -v x2=$length_x2 -v x3=$length_x3 'BEGIN{print sqrt(x1^2+x2^2+x3^2)}')
      length_y=$(awk -v y1=$length_y1 -v y2=$length_y2 -v y3=$length_y3 'BEGIN{print sqrt(y1^2+y2^2+y3^2)}')
      length_z=$(awk -v z1=$length_z1 -v z2=$length_z2 -v z3=$length_z3 'BEGIN{print sqrt(z1^2+z2^2+z3^2)}')
   #  remove old SPRING_R0 parameter
      sed -i '/SPRING_R0/d' ./INCAR
   #  add a new SPRING_R0 parameter   
      echo "SPRING_R0 = $length_x $length_y $length_z" >> INCAR
   fi 
#
#     For simulations with forces applied to one or several atoms 
#     in the system, update the SPRING_R0 parameters
#
   if [ $steered_opt == "-atom_steered" ]; then
      sed -nr '/SPRING_R0/p' REPORT > SPRING_lines
      tail -1 SPRING_lines > SPRING_last
      spring_all=`cat SPRING_last`   
      spring_words=( $spring_all )   
      spring_num=${#spring_words[@]} 

      sed -i '/SPRING_R0/d' ./INCAR
      echo "SPRING_R0 = ${spring_words[@]:2:spring_num}" >> INCAR
   fi	   


   #  Remove old ML_CTIFOR and ML_ISTART parameters, if present
   sed -i '/ML_CTIFOR/d' ./INCAR
   sed -i '/ML_ISTART/d' ./INCAR

   # Remove old TEBEG parameter 

   sed -i '/TEBEG/d' ./INCAR

   # Remove old NSW (step number) parameter
   sed -i '/NSW/d' ./INCAR

   # update the ML_CTIFOR parameter for the next trajectory part 
   ctifor_old=$ctifor
   ctifor=` grep "THRUPD" ML_LOGFILE | tail -1 | awk '{ printf( $4) }'`
   # If the variable is empty, use the previous value 
   if [ -z "$ctifor" ]; then
      ctifor=$ctifor_old
   fi	   
   echo "ML_CTIFOR = $ctifor" >> INCAR
   #echo "ML_ISTART = 1" >> INCAR

   # Add new starting temperature
   echo "TEBEG = $tebeg_act" >> INCAR

   # Add new MD step number
   echo "NSW = $nsw_act" >> INCAR
   
   sleep $sleept
done	

echo "  " >> ml_long.log
echo " The script 'ml_long.sh' has finished! Goodbye! " >> ml_long.log

