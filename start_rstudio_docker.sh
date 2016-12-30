#!/usr/bin/env bash

# https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image

# https://support.rstudio.com/hc/en-us/community/posts/202827628-Using-HTTPS-instead-of-HTTP
# https://www.digitalocean.com/community/tutorials/how-to-configure-nginx-with-ssl-as-a-reverse-proxy-for-jenkins

PASS="8_juggleD_albiNo_12_eleVens_cRush"
DOCKER_MNTPOINT="/home/rstudio/parker_rat_lung"
WORKSPACE_MNTPOINT="$DOCKER_MNTPOINT/workspace"
DATA_MNTPOINT="$DOCKER_MNTPOINT/raw_data"
HOST_BASE="$HOME/parker_rat_lung"
HOST_SCRATCH="/mnt/hts_scratch/Members/josh/parker_rat_lung"
WORKSPACE="$HOST_SCRATCH/workspace"
RAW_DATA="$HOST_SCRATCH/raw_data"
DOCKER_IMAGE_TAG="v5rc3"
DOCKER_IMAGE_NAME="granek/rstudio_qiime"
DOCKER_IMAGE="${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}"

PORT_NUMBER=8787
if [[ -z $3 ]] ; then
    echo "Using default PORT: $PORT_NUMBER"
else				
    PORT_NUMBER=$3
    echo "Using PORT from command line: '$PORT_NUMBER'"
fi

if [ "$1" == "shell" ]; then
    DOCKER_COMMAND="/bin/bash"
    CONTAINER_NAME="rstudio_shell_${DOCKER_IMAGE_TAG}"
    DOCKER_ARGS="--rm --interactive --tty --user rstudio"
    echo "------------------------------"
    echo "In docker run the following:"
    echo "cd $DOCKER_MNTPOINT"
    echo "make --dry-run fastqs -j4 -f analyze_culture_sequences.mk --warn-undefined-variables NUMTHREADS=4 RAW_DATA_DIR=$DOCKER_MNTPOINT/raw_data"
    echo "------------------------------"
elif [ "$1" == "root" ]; then
    DOCKER_COMMAND="/bin/bash"
    CONTAINER_NAME="rstudio_root_${DOCKER_IMAGE_TAG}"
    DOCKER_ARGS="--rm --interactive --tty"
elif [ "$1" == "rstudio" ]; then
    DOCKER_COMMAND=""
    CONTAINER_NAME="rstudio_web_${DOCKER_IMAGE_TAG}"
    DOCKER_ARGS="--detach --publish $PORT_NUMBER:8787 -e PASSWORD=$PASS"
else
   echo "Must supply command line argument! Should be run as one of the following commands:"
   echo "'$0 [shell|root|rstudio] (CONTAINER_NAME|default) (PORT_NUMBER)'"
   echo "The first argument (in []) is required, the second and third are optional"
   exit 1
fi

if [[ -z $2 ||  "$2" == "default" ]] ; then
    echo "Using default CONTAINER_NAME: $CONTAINER_NAME"
else				
    CONTAINER_NAME=$2
    echo "Using CONTAINER_NAME from command line: '$CONTAINER_NAME'"
fi



if [ ! -d "$WORKSPACE" ]; then
    echo "NOT MOUNTED: $WORKSPACE"
    echo "REFUSING TO START DOCKER"
    exit 1
fi

if [ ! -d "$RAW_DATA" ]; then
    echo "NOT MOUNTED: $RAW_DATA"
    echo "REFUSING TO START DOCKER"
    exit 1
fi

docker run $DOCKER_ARGS \
       --name $CONTAINER_NAME \
       -v $HOST_BASE:$DOCKER_MNTPOINT \
       -v $WORKSPACE:$WORKSPACE_MNTPOINT \
       -v $RAW_DATA:$DATA_MNTPOINT \
       $DOCKER_IMAGE \
       $DOCKER_COMMAND

# --------------------------------------------------
# GIT ENV variables (not sure if they are working)
# -e GIT_AUTHOR_NAME="Josh Granek" \
# -e GIT_AUTHOR_EMAIL="josh@duke.edu" \
# -e GIT_COMMITTER_NAME="Josh Granek" \
# -e GIT_COMMITTER_EMAIL="josh@duke.edu" \
# --------------------------------------------------
