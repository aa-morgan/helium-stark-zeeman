#!/bin/bash

# User variables
PYTHON_SCRIPT_NAME=legion_diagonalisation.py
PYTHON_SCRIPT_DIR=${HOME}/code/nlsjm/helium-stark-FS/notebooks
JOB_SCRIPT_NAME=serial_job_python.sh
JOB_SCRIPT_DIR=${HOME}/job_scripts/nlsjm
DATE_TIME=$(date +%F_%T | sed -e 's|:|-|g')
F_START=0.852
F_END=0.86
F_STEP=0.0002
B_FIELD=16.154E-4
N_MIN=72
N_MAX=73
S=1

# Create submitted temp folder for batch files
if [ ! -d $JOB_SCRIPT_DIR/submitted ]; then
mkdir $JOB_SCRIPT_DIR/submitted
fi
mkdir $JOB_SCRIPT_DIR/submitted/$DATE_TIME

# Submit job scripts
echo $(date) - Submitting jobs with fields, $F_START : $F_STEP : $F_END, and n from $N_MIN to $N_MAX
INDEX=1
for F in $(seq $F_START $F_STEP $F_END); do

# Make directory for single job
mkdir $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX

# Copy Python Script, and submission scripts into submitted temp folder
cp $PYTHON_SCRIPT_DIR/$PYTHON_SCRIPT_NAME $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX
cp $JOB_SCRIPT_DIR/batch_job.sh $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX
cp $JOB_SCRIPT_DIR/serial_job_python.sh $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX

# Edit python script with correct E field magnitude and vector
sed -ie "s|_Efield=.*|_Efield=${F}|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
sed -ie "s|_field_angle=.*|_field_angle=90.0|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
# Edit python script with correct B field
sed -ie "s|_Bfield=.*|_Bfield=${B_FIELD}|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
# Edit python script with correct nmin and nmax
sed -ie "s|_n_min=.*|_n_min=${N_MIN}|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
sed -ie "s|_n_max=.*|_n_max=${N_MAX}|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
# Edit python script with correct S
sed -ie "s|_S=.*|_S=${S}|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME
# Edit load/save directory
sed -ie "s|_matrices_dir=.*|_matrices_dir='data'|g" $JOB_SCRIPT_DIR/submitted/$DATE_TIME/$INDEX/$PYTHON_SCRIPT_NAME

# Submit the job
echo $(date) - Submit job with F $F, index $INDEX
# Parameters: N_MIN N_MAX DATE_TIME INDEX
qsub $JOB_SCRIPT_DIR/$JOB_SCRIPT_NAME $N_MIN $N_MAX $S $DATE_TIME $INDEX $PYTHON_SCRIPT_NAME

# Increment index
INDEX=$((INDEX+1))
done

echo $(date) - All done!
