#!/bin/sh

echo "manage_MDs.sh for VASP4CLINT by Andreas Mölkner, Version 3, Last Update: 30.04.2024"
echo ""
echo "This scirpt manages VASP ML_FF and AIMD MDs with different starting structures at different constant Temperatures!"
echo "The following Files have to be given in a Directory called input:"
echo "X POSCAR_files named POSCAR1, POSCAR2, ..."
echo "INCAR-file: TEBEG, TEEND and NSW will be adatped by the scirpt"
echo "POTCAR and KPOINTS"
echo "The MD.in file"
echo ""
echo "The script will be used as follows:"
echo "manage_MDs.sh [Action to be taken]"
echo ""
echo "Possible Actions:"
echo "0: Preparation of initial Directories"
echo "1: Preparation and Starting of calculations"
echo "2: Only Starting of the highest run-Number "
echo "3: Check if calculations finished"
echo "4: Generate new Directory for the continuation of the calcualtion (with NSW of ternary.in)"
echo "5: Restart unfinished calculations"
echo "6: Generate new Directory for the continuation of the calcualtion if the did not finish (with NSW beeing the difference)"
echo "7: Generate big combined XDATCAR of all runs for each Temperature" 

#This script has to be ran in the "MDs" 
echo ""
if [ -z "$1" ]; then
   echo "No command line argument submitted!"
   echo "Insert:"
   echo "0: Preparation of initial Directories"
   echo "1: Preparation and Starting of initial calculations"
   echo "2: Only Starting of the highest run-Number "
   echo "3: Check if calculations finished"
   echo "4: Generate new Directory for the continuation of the calcualtion (with NSW of ternary.in)"
   echo "5: Restart unfinished calculations"
   echo "6: Generate new Directory for the continuation of the calcualtion if the did not finish (with NSW beeing the difference)"
   echo "7: Generate big combined XDATCAR of all runs for each Temperature"
   exit 1
fi
action=$1

if [ "$action" != 0 ] && [ "$action" != 1 ] && [ "$action" != 2 ] && [ "$action" != 3 ] && [ "$action" != 4 ]  && [ "$action" != 5 ] && [ "$action" != 6 ]  && [ "$action" != 7 ] ;  then
   echo "Wrong command line argument! $action was submitted"
   echo "Insert:"
   echo "0: Preparation of Directories"
   echo "1: Preparation and Starting of calculations"
   echo "2: Only Starting"
   echo "3: Check for calculations"
   echo "4: Generate new Directory for continuation and start of calcualtion"
   echo "5: Restart unfinished calculations"
   echo "6: Generate new Directory for the continuation of the calcualtion if the did not finish (with NSW beeing the difference)"
   echo "7: Generate big combined XDATCAR of all runs for each Temperature"
   exit 2
elif [ "$action" == 0 ]; then
   echo "Only Preparation of Directories! Argument: $action" 
#   if [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
#      echo "Temperatures or Name of the System missing"
#      echo "manage_tern_MDs.sh [Action to be taken] [lowest Tempertaure] [highest Temperature] [dT] [Name of the System]"
#      exit 3
#   fi
elif [ "$action" == 1 ]; then
   echo "Preparation and Starting of Calcualtions! Argument: $action"
#   if [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
#      echo "Temperatures or Name of the System missing"
#      echo "manage_tern_MDs.sh [Action to be taken] [lowest Tempertaure] [highest Temperature] [dT] [Name of the System]"
#      exit 3
#   fi
elif [ "$action" == 2 ]; then
   echo "Only Starting of Calcualtions of the highest run-Number! Argument: $action"
elif [ "$action" == 3 ]; then
   echo "Check for calculations! Argument: $action"
elif [ "$action" == 4 ]; then
   echo "Generate new Directory for continuation and start of calcualtion! Argument: $action"
elif [ "$action" == 5 ]; then
   echo "Restart unfinished calcualtions! Argument: $action"
elif [ "$action" == 6 ]; then
   echo "Generate new Directory for the continuation of the calcualtion if the did not finish (with NSW beeing the difference)! Argument: $action"
elif [ "$action" == 7 ]; then
   echo "Generate big combined XDATCAR of all runs for each Temperature! Argument: $action"
fi

################### Read INPUT-file ################################

scalms=$( grep "name=" ./input/MD.in | awk '{print $2}' )
NTs=$( grep "NTs=" ./input/MD.in | awk '{print $2}' )
dt=$( grep "dt=" ./input/MD.in | awk '{print $2}' )

declare -a temp1
declare -a temp2
for ((i=1; i<=$NTs;i++)); do
   dummy=$( grep T$i= ./input/MD.in | awk  '{print $2}' )
   temp1+=($dummy)
   dummy2=$( grep T$i= ./input/MD.in | awk  '{print $3}' )
   temp2+=($dummy2)
done

NSW=$( grep "NSW=" ./input/MD.in | awk  '{print $2}' )

part=$( grep "partition=" ./input/MD.in | awk  '{print $2}' )
nodes=$( grep "nodes=" ./input/MD.in | awk  '{print $2}' )
cores=$( grep "cores=" ./input/MD.in | awk  '{print $2}' )
vasp=$( grep "VASP=" ./input/MD.in | awk  '{print $2}' )

echo ""
echo ""
echo "MD Settings:"
echo "Name= $scalms"
echo "NTs= $NTs"
echo "T_start: ${temp1[@]}"
echo "T_end:   ${temp2[@]}"
echo "NSW= $NSW"
echo "dt= $dt"

echo ""
echo "Slurm Settings:"
echo "Partition= $part"
echo "Nodes= $nodes"
echo "Cores= $cores"
echo "VASP= $vasp"
echo ""
echo ""
################# FUNCTIONS #####################################

function write_slurmscript {
rm slurm_script
touch slurm_script

echo "#! /bin/bash -l" >> slurm_script
echo "#" >> slurm_script 
echo "#SBATCH --nodes=$nodes" >> slurm_script
echo "#SBATCH --ntasks=$cores" >> slurm_script
echo "#SBATCH --time=24:00:00" >> slurm_script
echo "#SBATCH --job-name=$scalms-$1_$2" >> slurm_script
echo "#SBATCH --export=NONE" >> slurm_script
echo "#SBATCH --partition=$part" >> slurm_script
echo "#SBATCH --mail-user=andreas.moelkner@fau.de" >> slurm_script
echo "#SBATCH --mail-type=All" >> slurm_script 
echo "unset SLURM_EXPORT_ENV" >> slurm_script
echo "module load intel/2021.4.0 mkl/2021.4.0 intelmpi/2021.6.0 hdf5/1.10.7-impi-intel" >> slurm_script
if [ $vasp = "6.4.2" ]; then
	echo "VASP=/home/hpc/b146dc/b146dc11/Software/vasp6.4.2/vasp.6.4.2/bin/vasp_std" >> slurm_script
else
	echo "VASP=/home/hpc/b146dc/b146dc11/Software/vasp6.4.1/vasp.6.4.1/bin/vasp_std" >> slurm_script
fi
echo 'ldd $VASP' >> slurm_script
echo 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK' >> slurm_script
echo 'echo “OMP_NUM_THREADS=$OMP_NUM_THREADS”'  >> slurm_script
echo 'export OMP_PLACES=cores'  >> slurm_script
echo 'export OMP_PROC_BIND=true'  >> slurm_script
echo "touch CHGCAR" >> slurm_script
echo "touch WAVECAR"  >> slurm_script
echo 'srun $VASP > vasp.out 2> vasp.err'  >> slurm_script
echo 'rm CHG' >> slurm_script

}


function mod_incar {
sed -i '/TEBEG/d' ./INCAR
sed -i '/TEEND/d' ./INCAR
sed -i '/NSW/d' ./INCAR
sed -i '/NCORE/d' ./INCAR
sed -i '/POTIM/d' ./INCAR

echo "TEBEG = $1" >> INCAR
echo "TEEND = $2" >> INCAR
echo "NSW = $3" >> INCAR 
echo "NCORE = 26" >> ./INCAR
echo "POTIM = $dt" >> INCAR

if [ $part = "multinode" ]; then
	sed -i '/NCORE/d' ./INCAR
	echo "NCORE = 18" >> ./INCAR
fi

}


function count_dir {
items=$( ls -1 | wc -l )
declare -a runs
for ((n=1; n<=50;n++)); do
   if [ -d  $n  ] ; then
      runs+=($n)
   fi
done
echo "${runs[-1]}"
}

function dirlist { 
items=$( ls -1 | wc -l )
declare -a runs
for ((n=1; n<=50;n++)); do
   if [ -d  $n  ] ; then
      runs+=($n)
   fi
done
echo "${runs[@]}"
}




####################################################################


#Get the POSCAR Numbers for Loops and sort array

cd input
declare -a array
SUB='POSCAR'
for f in ./*; do
   if  [[ "$f" == *"$SUB"* ]];then
      NAME=${f##*CAR}
      array+=($NAME)
   fi
done
cd ../	
sorted=( $( printf "%s\n" "${array[@]}" | sort -n ) )



if [ "$action" == 0 ] || [ "$action" == 1 ]; then
   echo "Starting to prepare calculations" 
   echo "Preparing Calcualtions for $NTs Temperatures for each Structure:"

   for k in "${sorted[@]}"
   do
      mkdir $k
      for ((i=0; i<$NTs;i++))
      do
	 Tstart=${temp1[$i]}
	 Tend=${temp2[$i]}
	 x=$(( $i + 1 ))
	 mkdir ./$k/T$x 
	 cd ./$k/T$x
	 mkdir 1
	 cd ./1 
	 write_slurmscript $k T$x
	 cp ../../../input/POSCAR$k POSCAR
	 cp ../../../input/POTCAR .
	 cp ../../../input/KPOINTS .
	 cp ../../../input/INCAR .
	 cp ../../../input/ML_FF .
	 mod_incar $Tstart $Tend $NSW
	 cd ../../../ 
      done
   done
fi


if [ "$action" == 1 ] || [ "$action" == 2 ]; then
   echo "Starting calculation in highest numbered Directory!"
   echo ""

   for k in "${sorted[@]}"
   do
      cd $k
      for ((i=0; i<$NTs;i++))
      do
	 x=$(( $i + 1 ))	      
	 cd ./T$x
	 dir=$( count_dir )
	 cd $dir
	 echo "$k/T$x/$dir will be started!"
	 sbatch slurm_script
	 echo ""
         cd ../../
      done
      cd ../      
   done
fi


if [ "$action" == 3 ]; then
   echo "Checking for calcualtions in highest numbered Directory!"

   for k in "${sorted[@]}"
   do
      cd $k
      for ((i=0;i<$NTs;i++))
      do
	 x=$(( $i + 1 ))
         cd ./T$x
         dir=$( count_dir )
         cd $dir
	 NSWact=$( grep NSW INCAR | awk '{print $3}' )
	 step=$( grep T= vasp.out | tail -1 | awk '{print $1}' )
	 if [ $step -eq $NSWact ]; then
            echo "$k/T$x/$dir has finished! $step of $NSWact"
	 else
            echo "$k/T$x/$dir has NOT finished! $step of $NSWact"
	 fi
         cd ../../
      done
      cd ../
   done
fi


if [ "$action" == 4 ]; then
   echo "Creating new Directory for the next part of the Simulation"

   for k in "${sorted[@]}"
   do
      cd $k
      for ((i=0;i<$NTs;i++))
      do
         x=$(( $i + 1 ))
         cd ./T$x
         Tstart=${temp1[$i]}
         Tend=${temp2[$i]}
         dir=$( count_dir )
	 new_dir=$(($dir + 1))
	 echo "$new_dir"
	 mkdir $new_dir
	 cd $new_dir
	 cp ../$dir/INCAR .
         cp ../$dir/CONTCAR .
	 cp ../$dir/KPOINTS .
	 cp ../$dir/POTCAR .
	 cp ../$dir/ML_FF .
         write_slurmscript $k T$x
	 contcar2poscar
	 mod_incar $Tstart $Tend $NSW
	 cd ../../
      done
      cd ../
   done
fi


if [ "$action" == 5 ]; then
   echo "Checking for calcualtions in highest numbered Directory"
   echo "and restarting unfinished calculations in the same Directory"
   for k in "${sorted[@]}"
   do
      cd $k
      for ((i=0;i<$NTs;i++))
      do
         x=$(( $i + 1 ))
         cd ./T$x
         dir=$( count_dir )
         cd $dir
         NSWact=$( grep NSW INCAR | awk '{print $3}' )
         step=$( grep T= vasp.out | tail -1 | awk '{print $1}' )
         if [ $step -eq $NSWact ]; then
            echo "$k/T$x/$dir has finished! $step of $NSWact"
         else
            echo "$k/T$x/$dir has NOT finished! $step of $NSWact"
	    echo "Calcualtion will be restarted!"
	    sbatch slurm_script
         fi
         cd ../../
      done
      cd ../
   done
fi


if [ "$action" == 6 ]; then
   echo "Creating new Directory for the next part of the Simulation"

   for k in "${sorted[@]}"
   do
      cd $k
      for ((i=0;i<$NTs;i++))
      do
         x=$(( $i + 1 ))
         Tstart=${temp1[$i]}
         Tend=${temp2[$i]}
         cd ./T$x
         dir=$( count_dir )
	 cd $dir
         NSWact=$( grep NSW INCAR | awk '{print $3}' )
         step=$( grep T= vasp.out | tail -1 | awk '{print $1}' )
	 cd ../
         if [ $step -eq $NSWact ]; then
            echo "$k/T$x/$dir has finished! $step of $NSWact"
         else
            echo "$k/T$x/$dir has NOT finished! $step of $NSWact"
            echo "New Directory will be generatred:"
            new_dir=$(($dir + 1))
            echo "$new_dir"
            mkdir $new_dir
            cd $new_dir
            cp ../$dir/INCAR .
            cp ../$dir/CONTCAR .
            cp ../$dir/KPOINTS .
            cp ../$dir/POTCAR .
            cp ../$dir/ML_FF .
            write_slurmscript $k T$x
            contcar2poscar
            NSWact=$( grep NSW INCAR | awk '{print $3}' )
            step=$( grep T= ../$dir/vasp.out | tail -1 | awk '{print $1}' )
	    NSWnew=$(( $NSWact - $step ))
	    Tact=$(( $Tstart + ($Tend - $Tstart) * ( $NSW - $NSWnew )/$NSW ))
	    echo "Tact= $Tact"
            mod_incar $Tact $Tend $NSWnew
	    cd ../
	 fi
         cd ../
      done
      cd ../
   done
fi




if [ "$action" == 7 ]; then
   echo "Generation of combined XDATCAR"
   echo "This can and will take a while"
   echo "Happy Homenode"

   #Generate needed directories
   mkdir combined
   cd combined
   for ((i=0;i<$NTs;i++))
   do
      x=$(( $i + 1 ))
      mkdir T$x
   done
   cd ../
   

   #copy first XDATCAR with header
   for ((i=0;i<$NTs;i++))
   do
      x=$(( $i + 1 ))
      cd ${sorted[0]}
      cd T$x
      
      items=$( ls -1 | wc -l )
      declare -a runs
      for ((n=1; n<=50;n++)); do
      if [ -d  $n  ] ; then
         runs+=($n)
      fi
      done
    
      cd ${runs[0]}
      
      cp XDATCAR ../../../combined/T$x/

      cd ../../../
      unset runs
   done

   #append all other XDATCARS
   for ((i=0;i<$NTs;i++))
   do
      x=$(( $i + 1 ))
      for k in "${sorted[@]}"
      do
         cd $k/T$x
	 
         items=$( ls -1 | wc -l )
         declare -a runs
         for ((n=1; n<=50;n++)); do
         if [ -d  $n  ] ; then
         runs+=($n)
         fi
         done

         for d in "${runs[@]}"
	 do
            if [ $d -eq ${runs[0]} ] && [ $k -eq ${sorted[0]} ]  ; then
	       echo "$k/T$x/$d will be skiped"
	    else
	       cd $d
	       pwd
	       lines=$( wc XDATCAR | awk '{print $1}' )
               actlines=$(( $lines - 7 ))
               tail -$actlines XDATCAR  >> ../../../combined/T$x/XDATCAR
	       cd ../
	    fi
         done

	 cd ../../
	 unset runs
      done

   done


fi
