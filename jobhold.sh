#!/bin/bash
# Script credit : https://www.ccn.ucla.edu/wiki/index.php/How_to_have_a_script_wait_for_jobs_before_proceeding (Modified by me)
# jobhold.sh
# JB 8/2010
# usage: 
#       jobhold.sh "job_name <my_command>" (command and job name should have different cases or names)
# For commands that are submitted to cluster, hold until jobs have completed
#NOT TESTED FOR SLURM(sbatch)

shopt -s expand_aliases
declare -a full_cmd=($(echo "$@"))
jobname=${full_cmd[0]}
sleep_time=5 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
me=`whoami`

##CREATE tmp file, put cmd in tmp file and then execute the tmp file
TMP_DIR=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}') 
CLUSTER_OPTS=$(grep -i -w "cluster_options" parameters.txt | awk -F'=' '{print $2}')

tmp_file=$(mktemp -p $TMP_DIR)
touch $tmp_file
#printf -- "time ${full_cmd[@]/$1}" > $tmp_file
printf -- %s"\n"%s"\n" "#!/bin/bash" " $(echo ${full_cmd[@]/$jobname})" > $tmp_file
chmod a+x $tmp_file
#cat $tmp_file #debugging

CLUSTER="sh"
CLUSTER_QUEUE_CMD="(time sh $(echo $full_cmd | awk '{$1="";print}')) 1>> $TMP_DIR/$jobname.o 2>> $TMP_DIR/$jobname.e "
CLUSTER_JOBS_CMD="jobs -x"
if [[ $(builtin type -P "qsub" | wc -l) == 1 ]]; then 
    CLUSTER="SGE"
    CLUSTER_QUEUE_CMD="qsub $CLUSTER_OPTS -o $TMP_DIR/$jobname.o -e $TMP_DIR/$jobname.e -N $jobname $tmp_file" #time ${full_cmd[@]/$1}"
    CLUSTER_JOBS_CMD="qstat | grep $me"
elif [[ $(builtin type -P "sbatch" | wc -l) == 1 ]]; then 
    CLUSTER="SLURM"
    CLUSTER_QUEUE_CMD="sbatch $CLUSTER_OPTS --job-name $jobname $tmp_file" #time ${full_cmd[@]:1}"
    CLUSTER_JOBS_CMD="sstat | grep $me"
fi

#echo $CLUSTER
echo $CLUSTER_QUEUE_CMD
#echo $CLUSTER_JOBS_CMD

#sleep 5

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
echo "=Job $id($1)($tmp_file) started="
status=`myqstat | grep -o $id` # check to see if job is running
while [ -n "$status" ] # while $status is not empty
do
    echo "Job $id($1)($tmp_file)($TMP_DIR/$jobname) running..."
	status=`myqstat | grep $id`
    #status=`$CLUSTER_JOBS_CMD | grep $id`
    #echo $status
    #NEED to error handle here as well
    # if [[ $(echo $status | grep -q "Eqw") ]]; then
    #     echo "----Job $id($1)($tmp_file) ERROR---- "
    #     exit -2
    #     break;
    # fi
    sleep $sleep_time
done
echo "----Job $id($1)(sh:$tmp_file)(logs:$TMP_DIR/$jobname.[o/e]) complete---- "

cat $tmp_file >> $TMP_DIR/$jobname.o

rm $tmp_file

exit 0