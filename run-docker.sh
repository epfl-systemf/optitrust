#!/bin/bash
docker run -v $(pwd)/:/optitrust optitrust sh -c "cd /optitrust && tools/view_result.sh build_only $1 0"