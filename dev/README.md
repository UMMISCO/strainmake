# Using the Dockerfile

Build the Dockerfile into an image.

```sh
docker build -t <image_name> .
```

Create and run a container using the image. Add the git repo as a volume.

```sh
docker run --rm -it -v /path/to/repo:/project <image_name>
```

In the container, get the test data.

```sh
cd /project/.test
bash prepare.sh
```
