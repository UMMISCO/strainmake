To build and push the image, go into this folder (`.docker/`), then run the build script:

```sh
bash build.sh
```

The images are currently tagged based on `<date>.<time>`.
The images are built using the repository code at a given commit (see argument `PIPELINE_COMMIT_SHA` in [`Dockerfile`](Dockerfile)).