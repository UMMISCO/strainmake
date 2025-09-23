#!/usr/bin/env bash
# build and push strainmake Docker image

IMAGE=strainmake
TAG=$(date +%y%m%d.%H%M%S)

docker build -t bapt931894/${IMAGE}:latest -t bapt931894/${IMAGE}:${TAG} .
docker push bapt931894/${IMAGE}:latest
docker push bapt931894/${IMAGE}:${TAG}