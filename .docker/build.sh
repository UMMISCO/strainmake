#!/usr/bin/env bash
# build and push strainmake Docker image

IMAGE=strainmake
VERSION=v0.2.0-pre1

docker build -t bapt931894/${IMAGE}:latest -t bapt931894/${IMAGE}:${VERSION} .
docker push bapt931894/${IMAGE}:latest
docker push bapt931894/${IMAGE}:${VERSION}