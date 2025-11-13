# Installation

## Docker

### Run from an image

## Build from source

### Troubleshooting

#### Problem of libclang

libclang is needed by bindgen. There's a few ways to install it.

1. Use the system package manager (apt, dnf, pacman, etc.) to install `llvm` if accessible
2. Compile manually or automatically with `spack`. Here's an example for spack.

```
spack install cmake@3.27%gcc@4.8.5
spack install -j 4 llvm +clang ~lld ~lldb ~flang ~polly %gcc@13 ^cmake@3.27%gcc@4.8.5
```
After that, use `spack load llvm@18.1.3` to load the libclang, and `spack find -p llvm@18.1.3` to get the location.
Set `LIBCLANG_PATH` before proceeding the cargo build.
