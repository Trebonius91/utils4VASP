#!/bin/sh
#
#    prepare_slurm: Prepare efficient calculations of VASP
#      single points, being generated from vasp2trainset 
#      (md_traj=eval mode).
#      It copies the needed input (INCAR, KPOINTS , POTCAR)
#      from the main folder into the frame folders and 
#      writes a number of slurm_scripts, depending on the 
#      number of calculations that shall be done serially 
#      per slurm_script, in order to keep the number of 
#      single jobs in the queue manageable
#    Part of VASP4CLINT
#      Julien Steffen, 2025 (julien.steffen@fau.de)
#

VASP_PATH="/home/hpc/b146dc/b146dc10/programs/vasp6.4.1/vasp.6.4.1/bin/vasp_std"
VASP_MODULES="module load intel/2021.4.0 mkl/2021.4.0 intelmpi/2021.6.0 hdf5/1.10.7-impi-intel"
echo ""
echo "This script prepares the efficient calculation of many
   single points, after executing vasp2trainset -md_traj=eval.
   It copies the INCAR, KPOINTS and POTCAR files into the 
   single folders and generates slurm_script files that
   start a fraction of the jobs serially, such that the number
   of single jobs in the queue stays manageable. 

   Usage: prepare_slurm.sh [number of jobs per slurm_script]"
#
#    Check if all input files are present
#
if [ ! -e "POTCAR" ]; then
    echo "The POTCAR file is not there!"
    exit
fi

if [ ! -e "KPOINTS" ]; then
    echo "The KPOINTS file is not there!"
    exit
fi

if [ ! -e "INCAR" ]; then
    echo "The INCAR file is not there!"
    exit
fi

#
#    Obtain the number of jobs per slurm_script
#
if [ -n "$1" ]; then
   jobnum=$1
else
   echo " Please give the number of jobs as command line argument!"
   exit
fi

#
#    Count number of single point folders in current folder
#
while true; do
    folder="frame$((folder_num + 1))"
    if [ -d "$folder" ]; then
        folder_num=$((folder_num + 1))
    else
        break
    fi
done

echo " Number of frame folders: $folder_num"

#
#    Copy input files to all folders
#
echo " Copy VASP input files to folders ..."
for ((i=1; i<=folder_num; i++)); do
   cp POTCAR frame$i
   cp KPOINTS frame$i
   cp INCAR frame$i
   touch frame$i/slurm_script
#
#    Copy a single job slurm_script file to each folder, in case one 
#      of the job failed
#
    echo "#! /bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --partition=singlenode
#SBATCH --time=24:00:00
#SBATCH --job-name=singlepoint$i
#SBATCH --export=NONE

$VASP_MODULES
VASP=$VASP_PATH

mpirun -np 72 \$VASP > vasp.out 2> vasp.err
" > frame$i/slurm_script
done
echo " ... done!"

#
#    Calculate number of chunks
#
chunks=$(( (folder_num + jobnum - 1) / jobnum ))

echo " Number of chunks: $chunks"
echo " Number of folders per chunk: $jobnum"

#
#     Now generate slurm_scripts for each chunk, in which it is looped through
#     the folders it manages
#
for ((i=1; i<=chunks-1; i++)); do
    folder_first=$(( (i - 1) * jobnum ))
    echo "#! /bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --partition=singlenode
#SBATCH --time=24:00:00
#SBATCH --job-name=singlepoint_part$i
#SBATCH --export=NONE

$VASP_MODULES
VASP=$VASP_PATH

for ((i=1; i<=$jobnum; i++)); do
   cd frame\$(( $folder_first + i ))
   mpirun -np 72 \$VASP > vasp.out 2> vasp.err
   pwd
   cd ..
done
" > slurm_script$i
done

#
#    Slurm script for the remaining jobs, if modulo not zero
#
remainder=$(( folder_num - (chunks - 1) * jobnum  ))
echo " Number of folders in remaining chunk: $remainder"

folder_first=$(( (chunks - 1) * jobnum ))
    echo "#! /bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --partition=singlenode
#SBATCH --time=24:00:00
#SBATCH --job-name=singlepoint_part$chunks
#SBATCH --export=NONE

$VASP_MODULES
VASP=$VASP_PATH

for ((i=1; i<=$remainder; i++)); do
   cd frame\$(( $folder_first + i ))
   mpirun -np 72 \$VASP > vasp.out 2> vasp.err
   pwd
   cd ..
done
" > slurm_script$chunks


#
#    Generate a driver script to call all slurm_scripts
#
echo "#!/bin/sh
#
for ((i=1; i<=$chunks; i++)); do
   sbatch slurm_script\$i
   sleep 0.5
done
" > start_slurm.sh

echo " prepare_slurm.sh finished! Execute the start_slurm.script to"
echo "  start all slurm_scripts"
echo " If single calculations fail, they can still be restarted with"
echo "  single-job slurm_scripts in the respective folders."
echo " "
