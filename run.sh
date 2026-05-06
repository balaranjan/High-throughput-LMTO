#!/bin/bash

PROJECT_DIR="$(pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"
SRC_DIR="$PROJECT_DIR/src"
CHANGES_DIR="$PROJECT_DIR/changes/$(date +%Y%m%d_%H%M%S)"
VOLUME_NAME="${PROJECT_NAME}-workspace"
CONTAINER_NAME="${PROJECT_NAME}-dev"
WORKDIR="/home/duser/workdir"

mkdir -p "$CHANGES_DIR"

# Create volume if it doesn't exist, seed it with src/
if ! podman volume exists "$VOLUME_NAME"; then
  echo "🌱 [$PROJECT_NAME] Seeding workdir volume..."
  podman run --rm --pull=never \
    -v "$SRC_DIR":/seed:ro \
    -v "$VOLUME_NAME":"$WORKDIR" \
    "${PROJECT_NAME}:latest" \
    bash -c "cp -r /seed/. $WORKDIR/ && chown -R duser:duser $WORKDIR"
fi

# Start existing container if it exists, otherwise create new one
if podman container exists "$CONTAINER_NAME"; then
  echo "🔄 [$PROJECT_NAME] Reusing existing container..."
  # Remove stopped container and recreate to ensure /start.sh is used correctly
  podman rm "$CONTAINER_NAME"
  echo "🚀 [$PROJECT_NAME] Recreating container..."
  podman run -d --pull=never \
    --name "$CONTAINER_NAME" \
    -v "$VOLUME_NAME":"$WORKDIR" \
    -p 2222:22 \
    --shm-size=512m \
    "${PROJECT_NAME}:latest" \
    /start.sh
else
  echo "🚀 [$PROJECT_NAME] Creating new container..."
  podman run -d --pull=never \
    --name "$CONTAINER_NAME" \
    -v "$VOLUME_NAME":"$WORKDIR" \
    -p 2222:22 \
    --shm-size=512m \
    "${PROJECT_NAME}:latest" \
    /start.sh
fi

echo "✅ [$PROJECT_NAME] Container started."
echo "   SSH:    ssh -p 2222 duser@<host-ip>"
echo "   VSCode: Connect to ${PROJECT_NAME}-container (port 2222) → open $WORKDIR"
echo "   Stop:   podman stop $CONTAINER_NAME"
echo ""
echo "   Waiting for container to stop..."

# Block here until container stops
podman wait "$CONTAINER_NAME"

# On exit — diff volume contents vs original src
echo "🔍 [$PROJECT_NAME] Diffing changes..."

podman run --rm --pull=never \
  -v "$SRC_DIR":/original:ro \
  -v "$VOLUME_NAME":"$WORKDIR":ro \
  -v "$CHANGES_DIR":/changes \
  "${PROJECT_NAME}:latest" \
  bash -c "
    cd $WORKDIR
    find . -type f | while read f; do
      orig=\"/original/\$f\"
      if [ ! -f \"\$orig\" ]; then
        mkdir -p \"/changes/\$(dirname \$f)\"
        cp \"\$f\" \"/changes/\$f\"
        echo \"NEW: \$f\"
      elif ! diff -q \"\$f\" \"\$orig\" > /dev/null 2>&1; then
        mkdir -p \"/changes/\$(dirname \$f)\"
        cp \"\$f\" \"/changes/\$f\"
        echo \"MODIFIED: \$f\"
      fi
    done
  "

echo ""
echo "✅ [$PROJECT_NAME] Changed files saved to: $CHANGES_DIR"
echo "   Review and merge manually."
