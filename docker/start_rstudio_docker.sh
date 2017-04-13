#!/usr/bin/env bash

# https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image

# https://support.rstudio.com/hc/en-us/community/posts/202827628-Using-HTTPS-instead-of-HTTP
# https://www.digitalocean.com/community/tutorials/how-to-configure-nginx-with-ssl-as-a-reverse-proxy-for-jenkins

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DEFINITION_FILE="${DIR}/docker_setup.sh"
PASS_FILE="${DIR}/local_docker_setup.sh"
echo $DEFINITION_FILE
echo $PASS_FILE

source $DEFINITION_FILE

if [ -e  $PASS_FILE ] ; then
    source $PASS_FILE
    # $PASS_FILE should contain only the following definition, replacing "MYpassword" with a high quality password
    # PASS="MYpassword"

else
    PASS="albiNo_juggleD__cRu_12_eleVens_8_sh"
fi

echo $PASS


if [[ -z $3 ]] ; then
    echo "Using default PORT: $PORT_NUMBER"
else				
    PORT_NUMBER=$3
    echo "Using PORT from command line: '$PORT_NUMBER'"
fi

if [ "$1" == "shell" ]; then
    DOCKER_COMMAND="/bin/bash"
    CONTAINER_NAME="rstudio_shell_${DOCKER_IMAGE_TAG}"
    DOCKER_ARGS="--rm --interactive --tty -u $UID -v /etc/passwd:/etc/passwd:ro"
    echo "------------------------------"
    echo "In docker run the following:"
    echo "cd $DOCKER_MNTPOINT"
    echo "make --dry-run --warn-undefined-variables"
    # echo "make --dry-run fastqs -j4 -f analyze_culture_sequences.mk --warn-undefined-variables NUMTHREADS=4 RAW_DATA_DIR=$DOCKER_MNTPOINT/raw_data"
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

docker run $DOCKER_ARGS $TIME_ZONE \
       --name $CONTAINER_NAME \
       -e USERID=$UID \
       -v $HOST_BASE:$DOCKER_MNTPOINT \
       -v $WORKSPACE:$WORKSPACE_MNTPOINT \
       -v $RAW_DATA:$DATA_MNTPOINT \
       $DOCKER_IMAGE \
       $DOCKER_COMMAND

