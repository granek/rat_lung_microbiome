DOCKER_IMAGE_TAG="DIv1rc3"
DOCKER_IMAGE_NAME="granek/parker_rat_lung"
DOCKER_IMAGE="${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}"

DOCKER_MNTPOINT="$HOME/parker_rat_lung"
WORKSPACE_MNTPOINT="$DOCKER_MNTPOINT/workspace"
DATA_MNTPOINT="$DOCKER_MNTPOINT/raw_data"

PORT_NUMBER=8786

echo $DOCKER_IMAGE
