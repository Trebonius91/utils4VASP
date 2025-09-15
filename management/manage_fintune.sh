#!/bin/sh

echo "manage_finetune.sh by Andreas Mölkner"
echo ""
echo "This script manages the VASP singlepoint calculations for the fintuning procedure of a mace forcfield."
echo 'The following Files have to be given in a Directory called "input":'
echo '"INCAR", "KPOINTS", "POTCAR" and "slurmscript".'
echo ""
echo "The follwing Actions can be performed by setting the FIRST command line argument to:"
echo '0: Copy the input-files from "input" into the "frameX" directories.'
echo "1: Starting the calculations in the directories."
echo "2: Checking for finished/unfinished calculations"
echo ""
echo 'Groups of "frameX" directories on which the action will be applied can be selected by the SECOND command line argument:'
echo 'a: All "frameX" directories are selected (Default).'
echo 'b: A batch of the directories is selected via the thrid (start) and fourth (end) command line argument (e.g. 10 100, would select all frames from 10 to 100, 10 and 100 included).'
echo 's: A single frame is selected via the third command line argument.'
echo 'f: All "frameX" directories with a finished calculation are selected'
echo 'u: All "frameX" directories with a unfinished calculation are selected'
echo ""
echo ""
echo "#############################################################################################"
echo ""

if [ -z "$1" ]; then
   echo "No FIRST command line argument submitted!"
   echo "Insert:"
   echo '0: Copy the input-files from "input" into the "frameX" directories.'
   echo "1: Starting the calculations in the directories."
   echo "2: Checking for finished/unfinished calculations"
   exit 1
fi
action=$1

if [ "$action" != 0 ] && [ "$action" != 1 ] && [ "$action" != 2 ] ;  then
   echo "Wrong FIRST command line argument! $action was submitted!"
   echo "Insert:"
   echo '0: Copy the input-files from "input" into the "frameX" directories.'
   echo "1: Starting the calculations in the directories."
   echo "2: Checking for finished/unfinished calculations"
   exit 2
fi


if [ -z "$2" ]; then
   echo "No SECOND command line argument submitted!"
   echo 'All "frameX" directories will be selected!'
   selection=a
   echo "Small Delay to cancel the script"
   echo ""
   sleep 2
   echo ""
   echo "The script will continue!"
   echo ""
else
   selection=$2
fi

if [ "$selection" != a ] && [ "$selection" != b ] && [ "$selection" != s ] && [ "$selection" != f ] && [ "$selection" != u ] ;  then
   echo "Wrong SECOND command line argument! $selection was submitted!"
   echo "Insert:"
   echo 'a: All "frameX" directories are selected (Default).'
   echo 'b: A batch of the directories is selected via the thrid (start) and fourth (end) command line argument (e.g. 10 100, would select all frames from 10 to 100, 10 and 100 included).'
   echo 's: A single frame is selected via the third command line argument.'
   echo 'f: All "frameX" directories with a finished calculation are selected'
   echo 'u: All "frameX" directories with a unfinished calculation are selected'
   exit 3
fi

if [ "$selection" == b ] ; then
   if [ -z "$3" ] || [ -z "$4" ] ; then
      echo "No THIRD or FOURTH command line argument submitted!"
      exit 4
   else
      sta=$3
      end=$4
   fi
elif [ "$selection" == s ] ; then
   if [ -z "$3" ] ; then
      echo "No THIRD command line argument submitted!"
      exit 4
   else
      sel=$3
   fi
fi

echo "Action $action was selected!"

echo "Selection $selection was selected!"
echo ""
echo "#############################################################################################"
echo ""


#store the X from all frameX in a sorted array! ONLY THE NUMBERS!
declare -a array
SUB='frame'
for f in ./*; do
   if  [[ "$f" == *"$SUB"* ]];then
      NAME=${f##*me}
      array+=($NAME)
   fi
done
all=( $( printf "%s\n" "${array[@]}" | sort -n ) )

#screening the directories for finished and unfinished calucaltions


if [ "$selection" == u ] || [ "$selection" == f ] || [ "$action" == 2 ] ; then
   echo "Searching for finished and unfinished calculations!"
   finished=()
   unfinished=()
   : > finished
   : > unfinished
   for i in "${all[@]}"
   do
      file="./frame$i/vasp.out"
      if [[ -f $file ]] && grep -q "1 F=" "$file"; then
         echo "frame$i: Calculation finished" >> finished
	 finished+=("$i")
      else
         echo "frame$i: Calculation not finished" >> unfinished
	 unfinished+=("$i")
   fi
done

echo ""
echo "Number of finished:   ${#finished[@]}"
echo "Number of unfinished: ${#unfinished[@]}"
echo ""
   if [ "$action" == 2 ] ; then
      exit 5
   fi
fi

#Adjusting the List of frames to the selected ones


if [ "$selection" == a ] ; then
   list=("${all[@]}")
elif [ "$selection" == b ] ; then
   list=($(seq $sta $end))
elif [ "$selection" == s ] ; then   
   list=($sel)
elif [ "$selection" == f ] ; then
   list=("${finished[@]}")
elif [ "$selection" == u ] ; then
   list=("${unfinished[@]}")
fi



#MAIN LOOP over "list" 

for i in "${list[@]}"
do
   if [ "$action" == 0 ] ; then
      #echo "$i, action $action, selection $selection"
      cp input/* frame$i
   elif [ "$action" == 1 ] ; then
      #echo "$i, action $action, selection $selection"
      cd frame$i
      sbatch slurm_script
      cd ../
   fi
done
	
