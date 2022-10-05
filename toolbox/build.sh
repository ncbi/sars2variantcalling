#!/bin/bash -eu

progname=$0
pipeline=ncbi-sars2
version=v0.0.0
test=0

usage() {
cat <<EOF
--DESCRIPTION:
    Script to build a docker for a pipeline

--USAGE:
  $progname --version v<M.m.r>
    --version:    docker version tag, default ${version}
    --help:       to display this help
EOF
    exit 1
}

while (( $# > 0 ))
do
    case "x$1" in
        x--version  ) version="$2"  ; shift ;;
        x--test     ) test=1                ;;
        x--help     ) usage                 ;;
    esac
    shift
done


BRANCH=$(git rev-parse --abbrev-ref HEAD)
GITSHA=$(git rev-parse --short HEAD)
BUILDID=$(date '+%s')
DOCKERFILE=Dockerfile
PROJECT=$( curl -s "http://metadata.google.internal/computeMetadata/v1/project/project-id" -H "Metadata-Flavor: Google" )

docker image build -t $pipeline -f $DOCKERFILE .

if ((test)); then
    exit
fi

docker tag $pipeline us.gcr.io/${PROJECT}/$pipeline:$version
gcloud docker -a --verbosity=error
docker push us.gcr.io/${PROJECT}/$pipeline:$version

