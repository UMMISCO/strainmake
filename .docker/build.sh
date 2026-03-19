#!/usr/bin/env bash
# build and push strainmake Docker image

IMAGE=strainmake
VERSION=v0.2.0-pre1

docker build -t bapt931894/${IMAGE}:latest -t bapt931894/${IMAGE}:${VERSION} .

# verifying the image starts, reports a version, and lists commands
echo "==> Smoke-testing image..."
docker run --rm bapt931894/${IMAGE}:latest --version
docker run --rm bapt931894/${IMAGE}:latest --help

echo "==> Smoke-test passed. Pushing..."
docker push bapt931894/${IMAGE}:latest
docker push bapt931894/${IMAGE}:${VERSION}