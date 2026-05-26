# Installation

Currently, REST provides three installation methods: **Conda installation**, **using Docker or Singularity images**, and **building from source**. For beginner users (non-developers), we recommend Conda installation for its convenience and speed. For users with some programming experience or developers of programs and methods, building from source is recommended to achieve optimal computational performance and feature customization. For those who need to run REST on a variety of computing platforms, we recommend using Docker or Singularity images to ensure environment consistency and easy deployment.

## Conda Installation

REST has been successfully published on Anaconda, allowing users to install the REST executable (binary) directly via [Conda](https://docs.conda.org.cn/projects/conda/en/stable/user-guide/getting-started.html). Currently, we support Conda installation on Linux x86-64 and MacOS Arm64 systems. Windows users must use [WSL](https://learn.microsoft.com/zh-cn/windows/wsl/install) (Windows Subsystem for Linux) to install Conda and REST.

First, create a new conda environment (recommended):

```sh
conda create -n rest python=3.11 -c conda-forge
```

Here, `rest` is the name of the new environment (you may customize it), and `python=3.11` specifies Python version 3.11 (other Python versions are not currently supported). The `-c conda-forge` option specifies conda-forge as the installation channel. Then, activate the environment and install REST:

```sh
conda activate rest
conda install rest -c restgroup -c mokit -c conda-forge
```

The `-c` options specify the channels where the REST program and the external library `mokit` are hosted.

If you wish to install additional Python libraries or tools in the current environment (`rest`), use the `conda install [package_name] -c conda-forge` command. It is recommended to give priority to conda-forge as the installation channel to avoid dependency conflicts. We also recommend setting `conda-forge` as the default channel; see the [conda-forge documentation](https://forge.conda.org.cn/docs/user/introduction/) for details.

After installation, verify that REST is installed successfully (make sure the corresponding Conda environment is activated):

```sh
which rest 
rest -h
```

## Using Docker or Singularity Images

### Introduction to Container Technology

Container technology is an application-centric virtualization technology. Applications and all their dependencies are packaged into container images, which can then be copied and run on different hosts. Container technology effectively solves issues related to dependency library installation, software environment isolation, and environment portability, enabling rapid deployment.

Both [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity/) are container technologies that allow container images to run in exactly the same way on any system with the corresponding container engine installed, ensuring reproducibility, portability, and scalability of computational results across personal computers, local servers, large cloud platforms, and HPC clusters. Docker's strengths lie in its complete ecosystem and rich image resources, making it developer-friendly. Singularity is better suited for high-performance computing (HPC) scenarios, supports unprivileged user execution, and has better compatibility with cluster resource managers such as SLURM.

**Note**: The following Docker and Singularity commands assume that the corresponding programs are already installed on the system. Users should refer to the official documentation to install Docker or Singularity. Cluster users can contact their cluster administrators for technical support.

### Downloading Pre-built Docker/Singularity Images

We provide REST Docker/Singularity images for [online download](https://restgroup.github.io/igor_group/rest_download.html) on our group's [homepage](https://restgroup.github.io/igor_group), allowing users to obtain the latest image files (Docker Container Image / Apptainer (Singularity) Container Image).

- For downloaded Docker image archive files (`rest_[version].tar.gz`), load the provided REST image with:

  ```sh
  docker load -i rest_[version].tar.gz
  ```

  `version` denotes the image version, named as `year.minor` (e.g., `2025.01`, the current version). After loading, use `docker images` to verify the image was successfully loaded.

- Similarly, users can download Singularity image files (`rest_[version].sif`). The `.sif` file is itself a Singularity image and requires no additional loading. On systems with Singularity installed, users can directly use this file to run a Singularity container. Unless otherwise stated, Docker/Singularity image files with the same version number correspond to the same development progress of REST (Gitee main branch).

### Building REST Images

We provide Docker scripts (Dockerfile) for automatically creating REST images, available in our team's [Gitee repository](https://gitee.com/restgroup/rest_docker). Users can build REST images independently on a Docker-equipped system using these scripts (by cloning the repository via `git clone`). The command is:

```sh
docker build -t [name]:[version] -f Dockerfile .
```

Here, `name` and `version` are the image name and version you specify. `Dockerfile` is the corresponding Docker script file. In the current repository, the default `Dockerfile` points to `Dockerfile.abini`, which bundles REST with additional features such as the STM simulation program. For example, a build command using a specific script could be:

```sh
docker build -t rest:latest -f ./Dockerfile.abini .
```

After the command completes, use `docker images` to verify the image was successfully loaded.

### Converting Docker Images to Singularity

Users who already have REST Docker images can convert them to Singularity images locally and then upload them to servers or HPC clusters.

- **From an existing image**

  When Docker Daemon is available, run the following command to create a Singularity image from an existing Docker image:

  ```sh
  singularity build [name]_[version].sif docker-daemon://[name]:[version]
  ```

  The `name` and `version` after `//` should match the information shown by `docker images`. The `name` and `version` before them are the filename of the saved Singularity image — you may customize them, but matching the Docker image name and version is recommended for easier management.

- **From an image archive file**

  You can also create a Singularity image from a Docker image archive file (`.tar`):

  ```sh
  singularity build [name]_[version].sif docker-archive://[path_to_tar_file]
  ```

  **Note**: The default archive format provided by our team is `.tar.gz`. It is uncertain whether the above command directly supports this format; users may first decompress it to `.tar` format before conversion.

  If you have lost a previous image archive file, you can save an archive from a locally existing Docker image with:

  ```sh
  docker save -o [name]_[version].tar [name]:[version]
  ```

  You may customize the archive filename; the `name` and `version` after it should match the information shown by `docker images`, i.e., specifying the image to save.

## Building from Source

Users can pull the REST source code from the repository and, provided all external dependencies are satisfied, compile REST using Rust's package manager Cargo to obtain the executable. This approach maximizes REST's computational performance and allows feature customization through compile options. The team's [Gitee repository](https://gitee.com/restgroup/rest_workspace) provides an overview of the compilation workflow, including installation of the Rust compiler and Cargo, building external dependencies, and obtaining and compiling REST source code, for researchers interested in contributing to development.

Available build documentation and scripts for reference:
- Building with conda dependencies: see the [Developer Installation Guide](../contributor/compile-guide.md).
- Docker build workflow: see the [rest_docker](https://gitee.com/restgroup/rest_docker) project.
- Build-from-scratch workflow: see the `README.rst` file in [rest_workspace](https://gitee.com/restgroup/rest_workspace).
- Conda distribution build workflow: see the [rest-feedstock](https://github.com/RESTGroup/rest-feedstock) project.

<!-- ### Troubleshooting

#### Problem with libclang

libclang is required by `bindgen`. There are several ways to install it:

1. Use the system package manager (`apt`, `dnf`, `pacman`, etc.) to install `llvm` if accessible.
2. Compile manually or automatically with `spack`. Here's an example using spack:

```
spack install cmake@3.27%gcc@4.8.5
spack install -j 4 llvm +clang ~lld ~lldb ~flang ~polly %gcc@13 ^cmake@3.27%gcc@4.8.5
```

After that, use `spack load llvm@18.1.3` to load libclang, and `spack find -p llvm@18.1.3` to get the location. Set `LIBCLANG_PATH` before proceeding with the cargo build. -->