#!/bin/bash -l

# Batch script to run a serial job on Legion with the upgraded
# software stack under SGE.

# Parameters: N_MIN N_MAX S DATE_TIME INDEX PYTHON_SCRIPT

# User variables
N_MIN=$1
N_MAX=$2
S=$3
DATE_TIME=$4
INDEX=$5
ANGLE=90.0
PYTHON_SCRIPT_NAME=$6
PYTHON_SCRIPT_DIR=${HOME}/job_scripts/nlsjm/submitted/$DATE_TIME/$INDEX
IMPORT_FILE_DIR=${HOME}/data/nlsjm
IMPORT_FILE_MATRIX=n=${N_MIN}-${N_MAX}_L_max=None_S=${S}_MJ=None_MJ_max=None
IMPORT_FILE_S=stark_${IMPORT_FILE_MATRIX}_angle=${ANGLE}
IMPORT_FILE_Z=zeeman_${IMPORT_FILE_MATRIX}

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request X Hr:Min:Sec of wallclock time.
#$ -l h_rt=6:00:0

# 3. Request X gigabyte of RAM 
#$ -l mem=64G

# 4. Request X gigabyte of TMPDIR space (default is 10 GB) (Controlled on a per node basis)
#$ -l tmpfs=50G

# 5. Set the name of the job.
#$ -N hsfs

# 6. Specify that the nodes used should be exclusive to avoid other jobs causing a seg-fault
# -ac exclusive

# 6. Load modules
echo $(date) - Loading modules
module unload compilers
module load compilers/gnu/4.9.2
module load python3/recommended

# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
#$ -wd /home/ucapaam/Scratch/output

# 8. Copy the job and Python scripts, and precomputed data to $TMPDIR
echo $(date) - Copying $PYTHON_SCRIPT_NAME and submission script to $TMPDIR
cp $PYTHON_SCRIPT_DIR/$PYTHON_SCRIPT_NAME $TMPDIR
cp $PYTHON_SCRIPT_DIR/*.sh $TMPDIR

mkdir $TMPDIR/data
# Stark Data
echo $(date) - Copying ${IMPORT_FILE_S}.npz to $TMPDIR/data
cp $IMPORT_FILE_DIR/$IMPORT_FILE_S.npz $TMPDIR/data
# Zeeman Data
echo $(date) - Copying ${IMPORT_FILE_Z}.npz to $TMPDIR/data
cp $IMPORT_FILE_DIR/$IMPORT_FILE_Z.npz $TMPDIR/data

# 9. Your work *must* be done in $TMPDIR 
cd $TMPDIR

# 10. Run the application.
echo $(date) - Script $PYTHON_SCRIPT_NAME: Started 
python $PYTHON_SCRIPT_NAME
echo $(date) - Script $PYTHON_SCRIPT_NAME: Complete 

# 11. Create directory for job 
mkdir $HOME/Scratch/output/$JOB_ID

# Note: Can use Local2Scratch directive to automate the transfering of data from local, compute node TMPDIR to Scratch

# 12. Preferably, tar-up (archive) all output files onto the shared scratch area
# Delete large untar'd pre-computed file before moving to Scratch
#echo $(date) - Removing $TMPDIR/data/$IMPORT_FILE_S.npz from TMPDIR
#rm $TMPDIR/data/$IMPORT_FILE_S.npz
#echo $(date) - Removing $TMPDIR/data/$IMPORT_FILE_Z.npz from TMPDIR
#rm $TMPDIR/data/$IMPORT_FILE_Z.npz

echo $(date) - Compressing and copying $TMPDIR to $HOME/Scratch/output
tar zcvf $HOME/Scratch/output/$JOB_ID/files_from_job_$JOB_ID.tar.gz $TMPDIR

# 13. Move .e and .o files to job folder
echo $(date) - Moving job logs to job folder $JOB_ID
mv /home/ucapaam/Scratch/output/*.e$JOB_ID /home/ucapaam/Scratch/output/$JOB_ID/
mv /home/ucapaam/Scratch/output/*.o$JOB_ID /home/ucapaam/Scratch/output/$JOB_ID/

echo $(date) - All Done!
