# 安装

目前，REST 程序提供了3种方式进行安装，分别是 **Conda 安装**、**使用 Docker 或 Singularity 镜像**以及**从源码编译**。对于初级用户（非开发者），我们推荐使用 Conda 进行安装，方便快捷；对于有一定编程基础的用户或程序与方法开发者，建议从源码编译，以获得最佳的计算性能和功能定制；而对于需要在多种计算平台上运行 REST 的用户，我们推荐使用 Docker 或 Singularity 镜像，以实现环境的一致性和便捷的部署。

## 通过 Conda 安装

REST 程序已成功发布在 Anaconda，允许用户通过 [Conda](https://docs.conda.org.cn/projects/conda/en/stable/user-guide/getting-started.html) 直接安装 REST 可执行程序（二进制文件）。目前，我们支持 Linux x86-64 及 MacOS Arm64 系统下的 Conda 安装。对于 Windows 用户，推荐使用 [WSL2](https://learn.microsoft.com/zh-cn/windows/wsl/install) (Windows Subsystem for Linux) 安装 Conda 以及 REST 程序，以获得最稳定的使用体验。REST 也支持在 Windows 原生环境下通过 Conda 安装，但稳定性可能不如 Linux 版本。

具体的做法是，首先，创建一个新的 conda 环境（推荐）

```sh
conda create -n rest python=3.11 -c conda-forge
```

其中，`rest` 为新环境的名称，用户可自行指定，`python=3.11` 表示指定 Python 版本号为 3.11，目前暂不支持其他 python 版本。选项 `-c conda-forge` 指定了 conda-forge 作为安装管道（channel）。随后，激活该环境 (activate)，安装 REST 程序

```sh
conda activate rest
conda install rest -c restgroup -c mokit -c conda-forge
```

其中，`-c` 选项分别指定 REST 程序以及外部库 mokit 所在的安装管道（channel）。

若用户希望在当前环境（rest）中安装其他 python 库或工具，可通过 `conda install [package_name] -c conda-forge` 命令安装，推荐优先使用 conda-forge 作为安装管道（channel），以避免依赖冲突，推荐用户将 `conda-forge` 设置为默认的安装管道，具体做法可参考 [conda-forge 文档](https://forge.conda.org.cn/docs/user/introduction/)。

完成后，用户即可通过如下命令检查 REST 程序是否安装成功（确保运行时激活对应的 Conda 环境）

```sh
which rest 
rest -h
```

### Windows 用户说明

REST 程序支持在 Windows 系统上运行，有以下两种方式：

1. **WSL2（推荐）**：通过 WSL2 安装 Linux 环境后，按照上述 Linux 的 Conda 安装方式进行安装。WSL2 方式最为稳定，功能完整，是官方推荐的 Windows 使用方案。

2. **原生 Windows**：REST 的 Conda 包也提供了 Windows 原生版本，用户可在 Windows 命令行中直接使用 `conda install` 安装。但需要注意，Windows 原生版本的稳定性可能不如 Linux 版本。

对于需要自动化安装的用户，团队在 [rest_workspace](https://gitee.com/restgroup/rest_workspace) 项目的 `install_scripts` 目录下提供了 Windows 安装脚本 (`install-rest.bat` 和 `install-rest-wsl.bat`)，分别用于原生 Windows 安装和 WSL 环境下的安装。

## 使用 Docker 或 Singularity 镜像

### 容器技术简介

容器技术是一种以应用软件为中心的虚拟化技术。以应用软件为单元，将软件及所有的依赖打包成容器镜像，打包后的容器镜像可直接拷贝到不同的主机上运行。通过容器技术，可以很好地解决安装软件时，依赖库的安装问题、软件环境的隔离以及软件环境的移植问题，实现快速部署。

[Docker](https://www.docker.com/) 和 [Singularity](https://sylabs.io/singularity/) 都属于容器技术，允许容器镜像可以在任何安装了相应容器引擎的系统上，以完全相同的方式运行，确保从个人电脑、本地服务器到大型云端或超算平台的计算结果的可重现性、可移植性和可扩展性。相较而言，Docker 的优势在于完整的生态系统和丰富的镜像资源，对软件开发友好；而 Singularity 则更适合高性能计算（HPC）场景，支持无特权用户运行，并且与集群资源管理器（如 SLURM）兼容性更好。

**注意**：以下 docker 与 singularity 相关命令建立在系统中已安装有对应的程序。用户应参考相关程序的官方网站或文档，在操作系统中完成 Docker 或 Singularity 的安装，对于集群用户，可联系集群管理员咨询相关的技术支持。

### 下载已有 Docker/Singularity 镜像

我们在课题组的 [GitHub 主页](https://restgroup.github.io/igor_group) 上提供了 REST Docker/Singularity 镜像的[在线下载](https://restgroup.github.io/igor_group/rest_download.html)，供用户获取最新版本的镜像文件(Docker Container Image/Apptainer (Singularity) Container Image)。

- 对于下载得到的 Docker 镜像 (Docker Container Image) 压缩文件 (rest_[version].tar.gz)，用户可以通过如下命令载入（load）我们提供的REST镜像

  ```sh
  docker load -i rest_[version].tar.gz
  ```

  `version` 表示镜像的版本号，以年份加小版本号的方式命名，例如 `2025.01`，即当前版本。载入完成后，用户可以通过 `docker images` 命令查看镜像是否成功载入。

- 同样，用户也可以下载 Singularity 镜像 (Apptainer (Singularity) Container Image) 文件 (rest_[version].sif)，该文件本身即为 Singularity 镜像文件，无需额外载入，对于已经安装Singularity的操作系统下，用户可直接使用该文件运行 Singularity 容器。如无额外说明，版本号相同的 Docker/Singularity 镜像文件对应相同开发进度的 REST 程序（Gitee 主分支）。

### 构建 REST 镜像

我们提供了用于自动创建 REST 镜像的 Docker 脚本（Dockerfile），可见团队的 [Gitee 仓库](https://gitee.com/restgroup/rest_docker)。用户可以使用这些脚本（通过 git clone 仓库地址），在安装有 Docker 的系统上自主构建 REST 镜像。命令为

```sh
docker build -t [name]:[version] -f Dockerfile .
```

其中，`name` 和 `version` 表示镜像的名称和版本号，用户可根据需要自行指定。`Dockerfile` 为对应的 Docker 脚本文件。在目前的仓库中，默认的 `Dockerfile` 指向 `Dockerfile.abini`，集合了 REST 程序本体以及计算平台的其他一些功能，如 STM 模拟程序。因此，一个使用指定脚本的构建命令可以是

```sh
docker build -t rest:latest -f ./Dockerfile.abini .
```

在上述命令执行完毕后，同样通过 `docker images` 命令查看镜像是否成功载入。

### 由 Docker 转换 Sigularity 镜像

对于已有 REST Docker 镜像的用户，可以通过如下途径，在本地将 Docker 镜像转换为 Singularity 镜像，再上传至服务器或 HPC 集群使用。

- **已存在的镜像**

  当 Docker Daemon 可用时，运行如下命令，从已有的 Docker 镜像创建 Singularity 镜像

  ```sh
  singularity build [name]_[version].sif docker-daemon://[name]:[version]
  ```

  其中，`//` 后的 `name` 和 `version` 表示镜像的名称和版本号，应与 `docker images` 命令中显示的信息一致。前者的 `name` 和 `version` 为保存的 Singularity 镜像的文件名，原则上用户可自行指定，但建议与 Docker 镜像名称和版本号保持一致，方便管理和识别。

- **镜像压缩文件**

  用户也可以运行如下命令，从 Docker 镜像压缩文件（.tar）创建 Singularity 镜像

  ```sh
  singularity build [name]_[version].sif docker-archive://[path_to_tar_file]
  ```

  **注意**：团队默认提供的压缩文件格式为 .tar.gz，不确定上述命令是否直接适用于该格式，用户可先解压为 .tar 格式后再进行转换。

  若用户已丢失了先前的镜像压缩文件，可以通过如下命令，从本地已有的 Docker 镜像中保存获取一个镜像压缩文件：

  ```sh
  docker save -o [name]_[version].tar [name]:[version]
  ```

  同样，用户可以自行指定压缩文件的名称，而后者的 `name` 和 `version` 应与 `docker images` 命令中显示的信息一致，即指定需要保存的镜像。

## 从源码编译

用户可以从代码仓库中拉取 REST 程序源码，在确保 REST 程序依赖的外部库完整的前提条件下，使用 Rust 语言的包管理器 Cargo 编译 REST 源码，获得可执行文件。这样的好处在于可以最大程度的保证 REST 的计算性能，同时可以通过编译选项，实现功能化定制。目前，在团队 [Gitee 仓库](https://gitee.com/restgroup/rest_workspace) 下，提供了大致的程序编译流程，包括 Rust 编译器与包管理器 Cargo 的安装，REST 外部依赖的构建以及 REST 程序源码的获取与编译的完整流程，供有兴趣参与开发的研究者参考。

可供参考的编译文档或脚本：
- 结合 conda 依赖的源码编译：参考 [开发者安装指南](../contributor/compile-guide.md)。
- Docker 编译流程：参考 [rest_docker](https://gitee.com/restgroup/rest_docker) 项目。
- 从头编译流程：参考 [rest_workspace](https://gitee.com/restgroup/rest_workspace) 项目中的 README.rst 文件。
- conda 发行版的构建流程：参考 [rest-feedstock](https://github.com/RESTGroup/rest-feedstock) 项目。

```{note}
REST 程序的 MPI 并行功能在编译时可选。如果编译环境不支持 MPI（如 Windows 原生环境），可通过以下方式禁用 MPI 编译：
`cargo build --no-default-features -F dftd3,dftd4,geometric-pyo3`
禁用 MPI 后，程序仍可正常进行所有量子化学计算，但仅支持线程级并行（Rayon），不支持跨节点 MPI 并行。
```
