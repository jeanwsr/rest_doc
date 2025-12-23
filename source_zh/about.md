# 序

本文档旨在作为 **REST 软件**（Rust-based Electronic Structure Toolkit）的用户和开发者提供的综合指南。文档内容涵盖了 REST 的多个方面，包括安装说明、方法使用细节、用户可控参数，以及开发实践的相关指南。

我们假定读者具备一定的**量子化学基础**知识，并且基本熟悉 **Linux 操作系统**。对于开发者而言，建议您对 **Rust 编程语言**有一定的了解。无论你是希望快速上手 REST 的新用户，还是希望为项目做出贡献的有经验开发者，我们都希望本手册能够为你有效地理解、使用和开发 REST 提供必要的信息。

## REST 功能概览

REST 是一个使用 Rust 编程语言编写的开源量子化学软件包，旨在为分子体系提供高效且精确的电子结构计算，并支持多种先进的现代电子结构方法。REST 致力于成为计算化学领域中一个多功能、易用的研究工具。其主要特性包括：

- **先进的电子结构方法**
REST 实现了多种前沿电子结构方法，包括 Hartree–Fock（HF）、密度泛函理论（DFT）、二阶 Møller–Plesset 微扰理论（MP2）、随机相位近似（RPA），以及更具特色的 XYG3 型双杂化密度泛函（xDH） 和 重整化 xDH（R-xDH）。

- **高效的双电子积分计算**
REST 采用基于 Resolution of Identity (RI，也称为 Density Fitting, DF) 的高效算法，对高斯型轨道（GTO） 的双电子积分进行加速计算。

- **基于 Rust 的张量工具包**
REST 利用自主开发的基于 Rust 语言的张量运算工具包（参见 [rest_tensors](https://gitee.com/restgroup/rest_tensors) 和 [rstsr](https://gitee.com/restgroup/rstsr)，实现了对多维数组和高维张量运算的高效处理，这是高性能电子结构计算程序的关键。
此外，这些工具包也使得 REST 的模块化开发以及新方法的集成更加简便，从而鼓励社区贡献与功能扩展。

### REST 2025.01 版本新特性

- 关键词和默认设置的更改
  - 使用 `max_memory` 以 MB 为单位指定内存限制（张颖，祝震予）。
  - 添加 `start_mix_cycle` 控制混合初始猜测的生成方式（虞凌岳）。
  - 为 RI 计算添加默认辅助基组（张颖，王石嵘）。

- 新方法和功能
  - 实现用于从分子簇计算周期性体系能带结构的 RRS-PBC 方法（林子涵）。
  - 实现并丰富用于激发态计算的 GW-BSE 方法（高琪芮）。
  - 支持 DSD-DH 泛函以及自定义参数的 SCS-MP2 计算（颜文杰）。

- 程序性能改进和优化
  - 优化并支持 SCF 过程中 RI 积分的新算法（祝震予）。
  - 修复受限开壳层（Restricted Open shell）计算和 Yamaguchi 自旋投影中的错误并进行改进（虞凌岳）。
  - 对 rest_libcint 的 API 进行全面更新（祝震予）。

## 贡献与引用 REST

REST 程序是由复旦大学化学理论研究中心开发，在徐昕教授的领导下，由张颖教授担任首席开发者完成。

如果您在研究中使用了 REST，请引用以下出版物：

```text
Zhiyun Li, Tianyi Gao, Shirong Wang, Sheng Bi, Rulin Feng, Zhenyu Zhu, Yilin Zhao, Wenjie Yan, Lingyue Yu, Qirui Gao, Zihan Lin, Jianming Wu, Igor Ying Zhang, Xin Xu. REST: Embracing the Rust Programming Language for Modern Electronic Structure Theory. Chinese Journal of Chemical Physics. DOI: 10.1063/1674-0068/cjcp2510156
```

## 进一步阅读

您可以访问 [REST 项目主页](https://gitee.com/restgroup/rest) 或 [REST 研究小组网站](https://restgroup.github.io/igor_group) 获取有关 REST 的更多信息，包括最新新闻、更新和其他资源。
