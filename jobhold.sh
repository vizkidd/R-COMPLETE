#!/bin/bash
# jobhold.sh
# JB 8/2010
# usage: 
#       jobhold.sh "job_name <my_command>"
# For commands that are submitted to cluster, hold until jobs have completed
#NOT TESTED FOR SLURM(sbatch)
shopt -s expand_aliases
full_cmd=$@
sleep_time=300 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
me=`whoami`
CLUSTER="sh"
CLUSTER_QUEUE_CMD="time sh $(echo $full_cmd | awk '{$1="";print}')"
CLUSTER_JOBS_CMD="jobs -x"
if [[ $(builtin type -P "qsub" | wc -l) == 1 ]]; then 
    CLUSTER="SGE"
    CLUSTER_QUEUE_CMD="qsub -V -N $1 time ${full_cmd[@]:1}"
    CLUSTER_JOBS_CMD="qstat | grep $me"
elif [[ $(builtin type -P "sbatch" | wc -l) == 1 ]]; then 
    CLUSTER="SLURM"
    CLUSTER_QUEUE_CMD="sbatch --job-name $1 time ${full_cmd[@]:1}"
    CLUSTER_JOBS_CMD="sstat | grep $me"
fi

#echo $CLUSTER
#echo $CLUSTER_QUEUE_CMD
#echo $CLUSTER_JOBS_CMD

sleep 5

alias myqstat=$CLUSTER_JOBS_CMD
export myqstat
#stdout=`$@` # call the command and capture the stdout

if [ $CLUSTER == "sh" ]; then 
    nohup $CLUSTER_QUEUE_CMD  &> /dev/null & # call the command and capture the stdout
    id="$!"
elif [ $CLUSTER == "SGE" ]; then 
    stdout=`$CLUSTER_QUEUE_CMD` # call the command and capture the stdout
    #echo $stdout
    id=`echo $stdout | grep -o "[[:digit:]]" |  paste -sd '' ` #awk -F' ' '{print $3}'` # get the jobid
elif [ $CLUSTER == "SLURM" ]; then 
    stdout=`$CLUSTER_QUEUE_CMD` # call the command and capture the stdout
    id=`echo $stdout | grep -o "[[:digit:]]" |  paste -sd ''` #awk -F' ' '{print $3}'` # get the jobid #UNTESTED
fi
#echo $id
#echo $stdout 
echo "=Job $id started="
status=`myqstat | grep -o $id` # check to see if job is running
while [ -n "$status" ] # while $status is not empty
do
    sleep $sleep_time
    echo "Job $id running..."
	status=`myqstat | grep $id`
    #status=`$CLUSTER_JOBS_CMD | grep $id`
    #echo $status
done
echo "----Job $id complete---- "