#!/bin/bash
#
#SBATCH --job-name=pwmig_test_multi
#SBATCH --output=/scratch1/10208/chenboyin/frontera/out/pwmig_test_multi_%j.out
#SBATCH --error=/scratch1/10208/chenboyin/frontera/err/pwmig_test_multi_%j.err
#SBATCH --time=3:00:00
#SBATCH --nodes=3                   
#SBATCH --ntasks-per-node=1         
#SBATCH -p rtx

pwd
date
rm -rf $SCRATCH/mongodb/data

mkdir -p $SCRATCH/mongodb/data
mkdir -p $SCRATCH/mongodb/log

#SECTION 3:  Define some basic control variables for this shell
# this sets the working directory
# SCRATCH is an environment variable defined for all jobs on stempede2
WORK_DIR=$SCRATCH/pwmig_work
# This defines the path to the docker container file.
# like SCRATCH WORK2 is an environment variable defining a file system
# on stampede2
PWMIG_CONTAINER=$SCRATCH/containers/parallel_pwmig_latest.sif
MSPASS_CONTAINER=$SCRATCH/containers/mspass_latest.sif
# specify the location where user wants to store the data
# should be in either tmp or scratch
DB_PATH=$SCRATCH/mongodb
# the base for all hostname addresses
HOSTNAME_BASE='frontera.tacc.utexas.edu'
# Sets whether to use sharding or not (here sharding is turned on)
DB_SHARDING=false
# define database that enable sharding
SHARD_DATABASE="usarraytest"
# define (collection:shard_key) pairs
SHARD_COLLECTIONS=(
    "arrival:_id"
)
# This variable is used to simplify launching each container
# Arguments are added to this string to launch each instance of a
# container.  stampede2 uses a package called singularity to launch
# each container instances
APP_COM="apptainer run ${PWMIG_CONTAINER}"
# APP_COM="apptainer run ${MSPASS_CONTAINER}"


# Section 4:  Set up some necessary communication channels
# obtain the hostname of the node, and generate a random port number
NODE_HOSTNAME=`hostname -s`
echo "primary node $NODE_HOSTNAME"
LOGIN_PORT=`echo $NODE_HOSTNAME | perl -ne 'print (($2+1).$3.$1) if /c\d(\d\d)-(\d)(\d\d)/;'`
STATUS_PORT=`echo "$LOGIN_PORT + 1" | bc -l`
echo "got login node port $LOGIN_PORT"

# create reverse tunnel port to login nodes.  Make one tunnel for each login so the user can just
# connect to stampede2.tacc.utexas.edu
for i in `seq 4`; do
    ssh -q -f -g -N -R $LOGIN_PORT:$NODE_HOSTNAME:8888 login$i
    ssh -q -f -g -N -R $STATUS_PORT:$NODE_HOSTNAME:8787 login$i
done
echo "Created reverse ports on Frontera logins"


# Section 5:  Launch all the containers
# In this job we create a working directory on stampede2's scratch area
# Most workflows may omit the mkdir and just use cd to a working
# directory created and populated earlier
mkdir -p $WORK_DIR
cd $WORK_DIR

# start a distributed scheduler container in the primary node
APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
APPTAINERENV_MSPASS_ROLE=scheduler $APP_COM &

# get the all the hostnames of worker nodes
WORKER_LIST=`scontrol show hostname ${SLURM_NODELIST} | \
             awk -vORS=, -v hostvar="$NODE_HOSTNAME" '{ if ($0!=hostvar) print $0 }' | \
             sed 's/,$/\n/'`
echo $WORKER_LIST

# start worker container in each worker node
APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
APPTAINERENV_MSPASS_SCHEDULER_ADDRESS=$NODE_HOSTNAME \
APPTAINERENV_MSPASS_ROLE=worker \
mpiexec.hydra -n $((SLURM_NNODES-1)) -ppn 1 -hosts $WORKER_LIST $APP_COM &

if [ "$DB_SHARDING" = true ] ; then
    echo 'Using Sharding MongoDB'
    # extract the hostname of each worker node
    OLD_IFS=$IFS
    IFS=","
    WORKER_LIST_ARR=($WORKER_LIST)
    IFS=$OLD_IFS

    # control the interval between mongo instance and mongo shell execution
    SLEEP_TIME=15

    # start a dbmanager container in the primary node
    username=`chenboyin`
    for i in ${!WORKER_LIST_ARR[@]}; do
        SHARD_LIST[$i]="rs$i/${WORKER_LIST_ARR[$i]}.${HOSTNAME_BASE}:27017"
        SHARD_ADDRESS[$i]="$username@${WORKER_LIST_ARR[$i]}.${HOSTNAME_BASE}"
        SHARD_DB_PATH[$i]="$username@${WORKER_LIST_ARR[$i]}.${HOSTNAME_BASE}:/tmp/db/data_shard_$i"
        SHARD_LOGS_PATH[$i]="$username@${WORKER_LIST_ARR[$i]}.${HOSTNAME_BASE}:/tmp/logs/mongo_log_shard_$i"
    done

    APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
    APPTAINERENV_MSPASS_SHARD_DATABASE=${SHARD_DATABASE} \
    APPTAINERENV_MSPASS_SHARD_COLLECTIONS=${SHARD_COLLECTIONS[@]} \
    APPTAINERENV_MSPASS_SHARD_LIST=${SHARD_LIST[@]} \
    APPTAINERENV_MSPASS_SLEEP_TIME=$SLEEP_TIME \
    APPTAINERENV_MSPASS_ROLE=dbmanager $APP_COM &

    # ensure enough time for dbmanager to finish
    sleep 30

    # start a shard container in each worker node
    # mipexec could be cleaner while ssh would induce more complexity
    for i in ${!WORKER_LIST_ARR[@]}; do
        APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
        APPTAINERENV_MSPASS_SHARD_ID=$i \
        APPTAINERENV_MSPASS_DB_DIR=$DB_PATH \
        APPTAINERENV_MSPASS_SLEEP_TIME=$SLEEP_TIME \
        APPTAINERENV_MSPASS_CONFIG_SERVER_ADDR="configserver/${NODE_HOSTNAME}.${HOSTNAME_BASE}:27018" \
        APPTAINERENV_MSPASS_ROLE=shard \
        mpiexec.hydra -n 1 -ppn 1 -hosts ${WORKER_LIST_ARR[i]} $APP_COM &
    done

    # Launch the jupyter notebook frontend in the primary node.
    # Run in batch mode if the script was
    # submitted with a "-b notebook.ipynb"
    if [ $# -eq 0 ]; then
        APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
        APPTAINERENV_MSPASS_SCHEDULER_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_DB_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_DB_DIR=$DB_PATH \
        APPTAINERENV_MSPASS_SHARD_ADDRESS=${SHARD_ADDRESS[@]} \
        APPTAINERENV_MSPASS_SHARD_DB_PATH=${SHARD_DB_PATH[@]} \
        APPTAINERENV_MSPASS_SHARD_LOGS_PATH=${SHARD_LOGS_PATH[@]} \
        APPTAINERENV_MSPASS_DB_MODE="shard" \
        APPTAINERENV_MSPASS_ROLE=frontend $APP_COM
    else
        while getopts "b:" flag
        do
            case "${flag}" in
                b) notebook_file=${OPTARG};
            esac
        done
        APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
        APPTAINERENV_MSPASS_SCHEDULER_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_DB_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_DB_DIR=$DB_PATH \
        APPTAINERENV_MSPASS_SHARD_ADDRESS=${SHARD_ADDRESS[@]} \
        APPTAINERENV_MSPASS_SHARD_DB_PATH=${SHARD_DB_PATH[@]} \
        APPTAINERENV_MSPASS_SHARD_LOGS_PATH=${SHARD_LOGS_PATH[@]} \
        APPTAINERENV_MSPASS_DB_MODE="shard" \
        APPTAINERENV_MSPASS_ROLE=frontend $APP_COM --batch $notebook_file
    fi
else
    echo "Using Single node MongoDB"
    # start a db container in the primary node
    APPTAINERENV_MSPASS_DB_DIR=$DB_PATH \
    APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
    APPTAINERENV_MSPASS_ROLE=db $APP_COM &
    # APPTAINERENV_MSPASS_ROLE=db apptainer instance start --bind $SCRATCH/mongodb/data:/data/db $SCRATCH/containers/mongo_4.4.sif mongodb_instance && apptainer exec instance://mongodb_instance mongod --bind_ip_all --fork --logpath $SCRATCH/mongodb/log/mongod.log &
    
    # ensure enough time for db instance to finish
    sleep 10
  
    # Launch the jupyter notebook frontend in the primary node.
    # Run in batch mode if the script was
    # submitted with a "-b notebook.ipynb"
    if [ $# -eq 0 ]; then
    APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
    APPTAINERENV_MSPASS_SCHEDULER_ADDRESS=$NODE_HOSTNAME \
    APPTAINERENV_MSPASS_DB_ADDRESS=$NODE_HOSTNAME \
    APPTAINERENV_MSPASS_SLEEP_TIME=$SLEEP_TIME \
    APPTAINERENV_MSPASS_ROLE=frontend apptainer exec --bind $SCRATCH/pwmig_work:/test $PWMIG_CONTAINER bash -c "
        cd /test
        python test.py -cs $NODE_HOSTNAME
        python pwmig_testsuite_dataprep.py -cs $NODE_HOSTNAME
        python pwstack.py -cs $NODE_HOSTNAME
        python pwmig.py -cs $NODE_HOSTNAME
    "
        # APPTAINERENV_MSPASS_ROLE=frontend $APP_COM 
    else
        while getopts "b:" flag
        do
            case "${flag}" in
                b) notebook_file=${OPTARG};
            esac
        done
        APPTAINERENV_MSPASS_WORK_DIR=$WORK_DIR \
        APPTAINERENV_MSPASS_SCHEDULER_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_DB_ADDRESS=$NODE_HOSTNAME \
        APPTAINERENV_MSPASS_SLEEP_TIME=$SLEEP_TIME \
        APPTAINERENV_MSPASS_ROLE=frontend $APP_COM --batch $notebook_file
        
    fi
fi