#!/bin/bash

START_ID=5340939
END_ID=5341044

JOB_SCRIPT_DIR=$HOME/job_scripts/nlsjm

# Create temp folder for batch files
if [ ! -d $JOB_SCRIPT_DIR/extracted ]; then
echo ${JOB_SCRIPT_DIR}/extracted
mkdir ${JOB_SCRIPT_DIR}/extracted
fi

# Delete previous data and temp files 
rm -rf $JOB_SCRIPT_DIR/extracted/data
rm -rf $JOB_SCRIPT_DIR/extracted/temp

if [ ! -d $JOB_SCRIPT_DIR/extracted/temp ]; then
mkdir $JOB_SCRIPT_DIR/extracted/temp
fi
if [ ! -d $JOB_SCRIPT_DIR/extracted/data ]; then
mkdir $JOB_SCRIPT_DIR/extracted/data
fi

for ID in $(seq $START_ID 1 $END_ID); do

# Copy job folder if it exists
if [ -d $HOME/Scratch/output/$ID ]; then

echo $(date) - Uncompressing job ID $ID
# Uncompress job TMPDIR into extracted/temp folder
tar -C $JOB_SCRIPT_DIR/extracted/temp -zxvf $HOME/Scratch/output/$ID/files_from_job_$ID.tar.gz
echo $(date) - Copying StarkMapData for job ID $ID into data folder
# Copy the required stark map file from tmpdir
cp $JOB_SCRIPT_DIR/extracted/temp/tmpdir/job/$ID.undefined/data/starkMap* $JOB_SCRIPT_DIR/extracted/data

fi

done

