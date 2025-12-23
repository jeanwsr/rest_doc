# 快速开始

## 快速安装

如果您希望快速安装并运行 REST 程序，我们推荐使用以下方式之一：

- **Conda 安装**
  如果您使用的是 Linux 操作系统，或 WSL（Windows Subsystem for Linux），并且已经安装了 Conda 包管理器，可以通过以下命令安装 REST (确保系统 GLIBC >= 2.28)：

  ```sh
  conda create -n rest python=3.11 -c conda-forge 
  codna activate rest
  conda install rest -c restgroup -c mokit/label/cf -c conda-forge
  ```

  准备好 REST 程序输入卡（`ctrl.in`）后，您可以通过以下命令，在命令行运行：

  ```sh
  rest -i ctrl.in
  ```

- **Docker 镜像**
  若您获取了 REST 的 Docker 镜像压缩文件（rest_2025.01.tar.gz），通过以下命令加载镜像，并运行 REST （将输入卡置于当前目录 `pwd` 下）：

  ```sh
  docker load -i rest_2025.01.tar.gz
  docker run --rm -v $(pwd):/data rest:2025.01 bash -c "rest -i ctrl.in"
  ```

关于 REST 程序安装方式的更多细节，请参见[安装指南](./user/install.md)。

## 运行 REST 程序

### 命令行

在 Linux 或 MacOS 系统中，当编译完成，REST 程序的可执行文件生成在 `$REST_HOME/target/release/` 或 `$REST_HOME/target/debug/` 目录下，环境变量 `$REST_HOME` 表示 REST 程序主体所在的目录，通常为用户克隆（clone）[代码仓库](https://gitee.com/restgroup/rest_workspace) 时指定的目录，`target` 目录为 Cargo 默认的编译输出目录，`release` 和 `debug` 表示用户采用的编译模式。一个可能的路径是

```sh
$HOME/rest-workspace/target/release/rest
```

其中，`$HOME` 表示用户的主目录路径，最后的 `rest` 则为 REST 程序的可执行文件。用户可以将该路径添加至系统的环境变量 `PATH` 中，方便后续调用。

```sh
export PATH=$PATH:$HOME/rest-workspace/target/release
```

完成后，用户可以通过以下命令在命令行运行 REST 程序。

```sh
rest -i ctrl.in > job.out 
```

`ctrl.in` 为 REST 程序输入文件，通过 `-i` 指定输入文件的路径，当该选项缺省时，程序会在当前目录下寻找名为 `ctrl.in` 的文件，若不存在则报错退出。目前，REST 的标准输出为终端输出（屏幕打印），通过重定向符号 `>`，用户可以将标准输出保存至指定的文件中，例如 `job.out` 。

对于通过 Conda 安装的 REST 程序，`rest` 命令已被添加至 Conda 环境的 `PATH` 中。因此，您只需在激活对应的 Conda 环境后，直接运行 `rest` 命令即可，也可以通过 `which rest` 命令检查它的具体路径。

当用户需要批量运行 REST 进行计算时，可以编写 Shell 脚本，调用 REST 程序完成一系列的计算任务。一个简单的示例脚本 `run_rest_job.sh` 如下所示：

```sh
for inp in *.in; do
    # Extract the base name without extension
    base_name="${inp%.in}"
    # Run the REST job
    rest -i "$inp" > ${base_name}.log
done
```

该脚本遍历当前目录下的所有 `.in` 输入文件，调用 REST 计算并生成与输入文件同名的输出文件 `.log` 。

### Docker

- **命令行模式**

借助容器技术，通过 REST 的 Docker 镜像（见 [安装指南](./user/install.md) 的容器镜像部分），用户可在任何安装有 Docker 的系统上，以完全相同的方式运行 REST 程序。当我们启动一个 REST 的 Docker 容器时，实际上产生了一个可以执行 REST 程序的 linux 环境，当我们以 **命令行模式(Shell Mode)**  进入这个容器时，所需要做的和以上命令行运行 REST 程序的方式并没有什么区别。

以命令行模式运行镜像容器的命令为

```sh
docker run --rm -it [name]:[version] -w /opt /bin/bash 
```

选项 `--rm` 表示在运行结束时删除容器，若希望保留该容器（以方便下次更快速的启动），则不声明，同时可通过 `--name` 给容器命名；`-it` 表示两个选项，通常同时声明，表示以交互模式运行容器，在启动后为容器分配一个命令行; `-w` 选项用于声明进入容器后所在的容器内的路径，此处为 `/opt`; `/bin/bash` 代表启动后，在容器内执行 bash 命令。执行成功，预计显示为

```sh
admin@1a7dd53fd176:/opt$
```

表示进入容器后，用户名为 admin，位于 `/opt` 路径下，此时，我们可以直接使用 `rest` 命令运行程序（注意提供输入卡），因为在镜像打包时，它已被加入环境变量中。

当使用 Docker 运行 REST 时，用户的输入文件往往不存在于容器内部，为此，需要将外部的目录/文件映射至容器内，而后在容器内部执行相关的操作，退出容器后，就可以在原本的目录中得到运行后的结果（程序输出）。相应的命令为

```sh
docker run --rm -v /path/to/local/dir:/path/in/container -it [name]:[version] -w /opt /bin/bash 
```

通过 `-v` 选项，将存在于容器外部的本地路径 `/path/to/local/dir` 映射至容器内部的特定路径 `/path/in/container` 两者以 `:` 分隔，表示对应。例如，在最开始的示例中，我们将当前目录 `$(pwd)` 映射至容器内的 `/data` 路径下。这样，启动容器并切换至 `/data` 路径，就可以看到当前目录下的所有文件，包括输入卡 `ctrl.in`，从而运行 REST 程序。用户也可以同时映射多个本地目录至容器内的不同路径，此时，需以多个 `-v` 选项并列声明。

- **Exec 模式**

对于 Exec 模式，可以理解为进入容器后执行对应的命令，执行完毕后退出。例如以下命令:

```sh
docker run --rm -v $(pwd):/data \
    -it rest:2025.01 \ 
    -w /data \
    /bin/bash -c "rest -i ctrl.in"
```

相比于命令行模式，唯一的不同在于最后，命令 `bash -c` 表示从字符串中读取并执行命令。因此，完整命令的含义是，启动容器，进入 `/data`，并执行 `rest` 读入输入卡 `ctrl.in` 执行计算，输出到屏幕。对于一般的情况，对应的方式可以是

```sh
docker run --rm -v /path/to/local/dir:/path/in/container \
    -it [name]:[version] \ 
    -w /path/in/container \
    /bin/bash -c "rest -i test.in" > test.log
```

这样执行完毕，得到的是名为 `test.log` 的输出文件，我们假定输入文件 `test.in` 位于 `/path/to/local/dir`，这可以是相对路径也可以是绝对路径。而对于批量计算的情景，只需将最后一行替换为

```sh
/bin/bash run_rest_job.sh 
```

即执行一个批量计算的 shell 脚本即可，可以是我们在前文展示的单个目录下存有所有输入文件的场景 `run_set_job.sh`，或更复杂的存储结构，这需要对应的脚本逻辑。

### Singularity

Singularity 命令的逻辑与 Docker 是类似的，由于我们一般在服务器与集群上使用，这里我们主要介绍 `exec` 命令相关，它与 Docker `exec` 模式的选项相近，一个通常的命令为

```sh
singularity exec rest_[version].sif bash -c "rest"
```

这表示运行对应的 REST Singularity File `rest_[version].sif` 生成容器，进入容器后执行 `rest` 命令。与 Docker 不同，Singularity 默认会将一些本地目录映射入容器中，包括当前工作目录 `$pwd` 。若要额外映射其他目录，则通过 `--bind` 指定，比如将本地的基组文件映射至容器的 `/data` 下

```sh
singularity exec --bind ./basis_set_pool:/data ...
```
