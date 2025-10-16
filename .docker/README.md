To build and push the image, go into this folder (`.docker/`), then run the build script:

```sh
bash build.sh
```

The images are currently tagged based on the release version (e.g., `v0.1.0-pre1`).
The images are built using the repository code at a given version tag (see argument `PIPELINE_VERSION` in [`Dockerfile`](Dockerfile)).