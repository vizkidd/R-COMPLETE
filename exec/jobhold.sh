#!/bin/bash
# Script credit : https://www.ccn.ucla.edu/wiki/index.php/How_to_have_a_script_wait_for_jobs_before_proceeding (Modified by me)
# jobhold.sh
# JB 8/2010
# usage: 
#       jobhold.sh parameters_file "job_name <my_command>" (command and job name should have different cases or names)
# For commands that are submitted to cluster, hold until jobs have completed
#NOT TESTED FOR SLURM(sbatch)

## FUNCTIONS

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
        echo "** Trapped CTRL-C"
        paused_proc=($(ps -o pid,stat -C jobhold.sh | sed 1d | grep -i -w 't' | awk '{print $1}' | grep -v "[[:alpha:]]\|[[:punct:]]" | head -n $MAX_CONCURRENT_JOBS))
        if [[ ! -z $paused_proc ]]; then #$CLUSTER == "sh" &&
            #Resume few stopped jobs
            echo ${paused_proc[@]}
            if [[ $(which parallel | awk '{print NR}') != 0 ]]; then
                parallel "kill -CONT {}" ::: $(echo ${paused_proc[@]: 0:$MAX_CONCURRENT_JOBS}) #${paused_proc[@]} 
            else
                for p_proc in $(echo ${paused_proc[@]: 0:$MAX_CONCURRENT_JOBS}) ; do
                    kill -CONT $p_proc
                done
            fi
            
        fi
        exit 255
}

function check_param() {
  local STDIN=$(cat)
  local tmp_name=$(echo $STDIN | awk -F'==' '{print $1}')
  local tmp_var=$(echo $STDIN | awk -F'==' '{print $2}')
  if [[ -z $tmp_var ]]; then
    echo "Parameter Error : $tmp_name is empty"
    exit 1
  fi
  echo $tmp_var
}

##ENTRYPOINT

shopt -s expand_aliases
declare -a full_cmd=($(echo "$@"))
#declare -a full_cmd=($@)
jobname=${full_cmd[1]}
only_cmd=$(echo "${full_cmd[@]: 2:${#full_cmd[@]}}")
sleep_time=10 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
me=`whoami`

if [[ -z ${full_cmd[0]} || ! -s ${full_cmd[0]} ]] ; then
  #echo "ERROR: Missing ${full_cmd[0]}!"
  #echo "Need : ${full_cmd[0]}"
  #exit -1
  if [[ ! -z $TMPDIR ]]; then
    TMP_DIR=$(echo $TMPDIR)
  else
    TMP_DIR="/tmp/"
  fi
  CLUSTER_OPTS=""
  MAX_CONCURRENT_JOBS=1
else
    TMP_DIR=$(grep -i -w "temp_path" ${full_cmd[0]} | check_param) 
    CLUSTER_OPTS=$(grep -i -w "cluster_options" ${full_cmd[0]} | check_param)
    MAX_CONCURRENT_JOBS=$(grep -i -w "max_concurrent_jobs" ${full_cmd[0]} | check_param)
fi

##CREATE tmp file, put cmd in tmp file and then execute the tmp file

# if [[ -z $TMP_DIR || -z $CLUSTER_OPTS || -z $MAX_CONCURRENT_JOBS ]] ; then
#   echo "ERROR: Missing parameters! (temp_path, cluster_options, max_concurrent_jobs)"
#   echo "Check ${full_cmd[0]}"
#   exit -1
# fi

mkdir -p $TMP_DIR
tmp_file=$(mktemp -p $TMP_DIR)
touch $tmp_file
#printf -- "time ${full_cmd[@]/$1}" > $tmp_file
#printf -- %s"\n"%s"\n" "#!$SHELL" " $(echo ${full_cmd[@]/$jobname})" > $tmp_file
printf -- %s"\n"%s"\n" "#!$SHELL" "$only_cmd" > $tmp_file
#printf -- %s"\n"%s"\n" "#!/bin/bash" $only_cmd > $tmp_file
chmod a+x $tmp_file
#cat $tmp_file #debugging

#$$ is current shell's process ID, we can get the parent shell of current shell and see how many child processes it has
#CURRENT_PARENT_SHELL=(ps -fC jobhold.sh) #$(ps -ho ppid --pid $$ | awk '{print $1}') 
CLUSTER="$SHELL"
CLUSTER_QUEUE_CMD="$SHELL $tmp_file" #CLUSTER_QUEUE_CMD="(time sh $(echo $tmp_file | awk '{$1="";print}')) 1>> $TMP_DIR/$jobname.o 2>> $TMP_DIR/$jobname.e "
CLUSTER_JOBS_CMD="jobs -rl" #CLUSTER_JOBS_CMD="jobs -x"
PROCS_COUNT_CMD="ps -o pid,stat -C jobhold.sh" #PROCS_COUNT_CMD="ps -ho pid,stat --ppid $CURRENT_PARENT_SHELL"
#PROCS_COUNT="jobs -rl"  
if [[ $(builtin type -P "qsub" | wc -l) == 1 ]]; then 
    CLUSTER="SGE"
    CLUSTER_QUEUE_CMD="qsub $CLUSTER_OPTS -o $TMP_DIR/$jobname.o -e $TMP_DIR/$jobname.e -N $jobname $tmp_file" #time ${full_cmd[@]/$1}"
    CLUSTER_JOBS_CMD="qstat -u $me" #CLUSTER_JOBS_CMD="qstat | grep $me"
    PROCS_COUNT_CMD="qstat -u $me" # | awk 'END{print NR-2;}'
elif [[ $(builtin type -P "sbatch" | wc -l) == 1 ]]; then 
    CLUSTER="SLURM"
    CLUSTER_QUEUE_CMD="sbatch $CLUSTER_OPTS --job-name $jobname $tmp_file" #time ${full_cmd[@]:1}"
    CLUSTER_JOBS_CMD="sstat -u $me" #CLUSTER_JOBS_CMD="sstat | grep $me"
    PROCS_COUNT_CMD="sstat -u $me" #| awk 'END{print NR-2;}'
fi

#sleep $[ ( $RANDOM % 10 )  + 1 ]s

PROCS_COUNT=($($PROCS_COUNT_CMD | sed 1d | grep -i -w 'r\|s' | awk '{print $1}'| grep -v "[[:alpha:]]\|[[:punct:]]")) #($(echo $($PROCS_COUNT_CMD) | awk '{print $1}'))
JOBHOLD_COUNT=${#PROCS_COUNT[@]}
#echo $CURRENT_PARENT_SHELL
#echo $PROCS_COUNT_CMD
#echo $JOBHOLD_COUNT
#echo ${PROCS_COUNT[@]}
#ps -o pid,stat,cmd -C jobhold.sh

while [[ "$JOBHOLD_COUNT" -gt $(echo "$MAX_CONCURRENT_JOBS"+1 | bc) ]];
do
  printf %s"\n" "Waiting for jobs to end..($JOBHOLD_COUNT)"
  kill -STOP $$
  #sleep $sleep_time
#  if [ $CLUSTER == "$SHELL" ]; then 
#    kill -STOP $$
#elif [ $CLUSTER == "SGE" ]; then 
#    sleep $sleep_time
#elif [ $CLUSTER == "SLURM" ]; then 
#    sleep $sleep_time
#fi

  PROCS_COUNT=($($PROCS_COUNT_CMD | sed 1d | grep -i -w 'r\|s' | awk '{print $1}'| grep -v "[[:alpha:]]\|[[:punct:]]"))
  JOBHOLD_COUNT=${#PROCS_COUNT[@]}
done

#echo $CLUSTER
echo $CLUSTER_QUEUE_CMD
#echo $CLUSTER_JOBS_CMD

#sleep 5

alias myqstat=$CLUSTER_JOBS_CMD
export myqstat
#stdout=`$@` # call the command and capture the stdout

if [ $CLUSTER == "$SHELL" ]; then 
    nohup $CLUSTER_QUEUE_CMD 1>> $TMP_DIR/$jobname.o 2>> $TMP_DIR/$jobname.e& #&> /dev/null & # call the command and capture the stdout
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

echo "=Job $id($jobname)($tmp_file) started="
#echo $(jobs -rl)
status=`myqstat | grep -o $id` # check to see if job is running
while [ -n "$status" ] # while $status is not empty
do
    PROCS_COUNT=($($PROCS_COUNT_CMD | sed 1d | grep -i -w 'r\|s' | awk '{print $1}' | grep -v "[[:alpha:]]\|[[:punct:]]")) #PROCS_COUNT=($($PROCS_COUNT_CMD | grep -i 'r' | awk '{print $1}'))
    JOBHOLD_COUNT=${#PROCS_COUNT[@]}
    #echo "Job $id($1)($tmp_file)($TMP_DIR/$jobname)($JOBHOLD_COUNT) running..."
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

printf %s"\n" "----Job $id($jobname)($SHELL:$tmp_file)(logs:$TMP_DIR/$jobname.[o/e]) complete---- "

#declare -x -g -i JOBHOLD_COUNT=$(($JOBHOLD_COUNT-1)) #export -n JOBHOLD_COUNT=$(($JOBHOLD_COUNT-1))

cat $tmp_file >> $TMP_DIR/$jobname.o

rm $tmp_file

#echo $($PROCS_COUNT_CMD)
paused_proc=($(ps -o pid,stat -C jobhold.sh | sed 1d | grep -i -w 't' | awk '{print $1}' | grep -v "[[:alpha:]]\|[[:punct:]]" | head -n $MAX_CONCURRENT_JOBS))
if [[ ! -z $paused_proc ]]; then #$CLUSTER == "sh" &&
    #Resume few stopped jobs
    echo ${paused_proc[@]}
    if [[ $(which parallel | awk '{print NR}') != 0 ]]; then
        parallel "kill -CONT {}" ::: $(echo ${paused_proc[@]: 0:$MAX_CONCURRENT_JOBS}) #${paused_proc[@]} 
    else
        for p_proc in $(echo ${paused_proc[@]: 0:$MAX_CONCURRENT_JOBS}) ; do
            kill -CONT $p_proc
        done
    fi
    
fi

exit 0