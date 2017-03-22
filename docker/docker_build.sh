#!/usr/bin/env bash                                                                                                                                                                      
# Find the directory of this script http://stackoverflow.com/a/246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source ${DIR}/docker_setup.sh
docker build -t ${DOCKER_IMAGE} $DIR

