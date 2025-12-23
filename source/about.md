# About

This manual is written as a comprehensive guide for both users and developers of the REST software, i.e., Rust based Electronic Structure Toolkit. This document aims to conver various aspects of REST, including installation instructions, usage details on methods and user-controlled parameters as well as guidelines for development practices. It is assumed that the reader can have a basic konwledge of Quantum Chemistry and is somehow familiar with Linux. For developers, some familiarity with the Rust programming language is recommended. Whether you are a new user looking to get started with REST or an experienced developer seeking to contribute to the project, we hope this manual will provide the necessary information to help you navigate and utilize the software effectively.

## REST Features Overview

REST is an open-source quantum chemistry software package written in Rust programming language. It is designed to provide efficient and accurate electronic structure calculations for molecular systems with advanced modern electronic structure methods. REST aims to be a versatile and user-friendly tool for researchers in the field of computational chemistry. Some of the key features of REST include:

- **Advanced Electronic Structure Methods**: REST implements a variety of state-of-the-art electronic structure methods, including Hartree-Fock (HF), Density Functional Theory (DFT), Second Order Møller-Plesset perturbation theory (MP2), Random Phase Approximation (RPA) and more distinctively, the XYG3 type of Doubly Hybrid Density Functionals (xDHs) as well as the Renormalized xDHs (R-xDHs).

- **Efficient Electronic Integrals**: REST utilizes efficient Resolution of Identity (RI, also known as Density Fitting, DF) based algorithms for the computation of electronic integrals with Gaussian Type Orbitals (GTO) to accelerate calculations.

- **Rust based Tensor Toolkit**: REST leverages custom-built tensor toolkits written in Rust (see [rest_tensors](https://gitee.com/restgroup/rest_tensors) and [rstsr](https://gitee.com/restgroup/rstsr)), providing efficient handling of multi-dimensional arrays and tensor operations, which is crucial for high performance electronic structure calculations. Moreover, with these toolkits at hand, modular development and integration of new methods with REST becomes more straightforward, encouraging community contributions and extensions.

### Features in REST 2025.01

- Changes to keywords and default behaviors
  - Use `max_memory` to specify memory limit in MB (Igor Ying Zhang; Zhenyu Zhu).
  - Add `start_mix_cycle` to control how mixed initial guess is generated (Lingyue Yu).
  - Add default auxiliary basis set for RI calculations (Igor Ying Zhang; Shirong Wang).

- New methods and functionalities
  - Implement RRS-PBC method for calculating band structures of periodic systems from molecular clusters (Zihan Lin).
  - Implement and enrich the GW-BSE method for excited state calculations (Qirui Gao).
  - Support for DSD-DH functionals and user-defined parameters in spin-component-scaled (SCS) MP2 calcualtions (Wenjie Yan).

- Performance improvements and optimizations
  - Optimize and support new algorithms of RI integrals during SCF (Zhenyu Zhu).
  - Bug fix and improvement in restricted open-shell calculations and Yamaguchi's spin projection (Lingyue Yu).
  - General update in APIs to rest_libcint (Zhenyu Zhu).

## Citing REST

If you use REST in your research, please cite the following publication:

```text
Zhiyun Li, Tianyi Gao, Shirong Wang, Sheng Bi, Rulin Feng, Zhenyu Zhu, Yilin Zhao, Wenjie Yan, Lingyue Yu, Qirui Gao, Zihan Lin, Jianming Wu, Igor Ying Zhang, Xin Xu. REST: Embracing the Rust Programming Language for Modern Electronic Structure Theory. Chinese Journal of Chemical Physics. DOI: 10.1063/1674-0068/cjcp2510156
```

## Further readings

You can visit either the [REST project page](https://gitee.com/restgroup/rest) or the [REST group website](https://restgroup.github.io/igor_group) for more information about REST, including the latest news, updates, and additional resources.
