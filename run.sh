#!/bin/sh
# Run a command inside the pi4decoder-dev container with the repo (and the
# test data directories) mounted. Usage: ./run.sh make test
exec docker run --rm --user "$(id -u):$(id -g)" --name "pi4decoder-${PI4_TASK:-run}-$$" --label owner=claude \
  -v /data/claude/pi4decoder:/work -w /work ${PI4_EXTRA_MOUNTS} pi4decoder-dev "$@"
