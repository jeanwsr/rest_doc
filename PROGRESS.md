# REST 文档更新进度

基于近期合并的 PR 扫描 (`gitee pr list --repo restgroup/rest --state merged`)，记录文档更新状态。

扫描日期：2026-06-18

## 已完成

| 项目 | 相关 PR | 文件 | 说明 |
|---|---|---|---|
| Smearing 文档 | #128 | `source_zh/user/scf.md` | 已由 !9 合入 rest_doc |
| PCM/溶剂化文档 | #136, #125 | `source_zh/user/solvent.md` | 新建，涵盖 4 种 PCM 模型、5 个关键词、空腔构建、示例 |
| 构型优化文档 | #134 | `source_zh/user/geomopt.md` | 重写，涵盖 geomeTRIC/L-BFGS 引擎、固定原子、TS/IRC、全部参数 |
| Windows 支持 | #126 | `source_zh/user/install.md` | 新增 Windows 用户说明（WSL2 推荐 + 原生安装） |
| MPI 可选编译 | #122 | `source_zh/user/install.md` | 新增 MPI 可选编译说明 |
| 溶剂化索引 | — | `source_zh/user/index.md` | toctree 中新增 `solvent` 条目 |

## 待更新

| 项目 | 相关 PR | 优先级 | 说明 |
|---|---|---|---|
| EDIIS/ADIIS mixer 文档 | #131 | 中 | `scf.md` 中 mixer 部分仅有 direct/linear/diis/ddiis；EDIIS 已写入 PR !150 README，待 !150 合并后同步至 `scf.md` |
| 外电场及解析梯度 | #145, #146 | 中 | `ext_field_dipole` 关键词未在用户文档中说明（仅在 README `[geom]` 区块），且梯度支持为新功能 |
| 118 元素 + IUPAC 2021 原子量 | #133 | 低 | 可在 `about.md` 或版本更新记录中简要提及 |
| VXC 格点优化关键词 | #127 | 低 | `vxc_screen_threshold`、`ao_cutoff`、`non0tab_blksize`、`drop_dense_ao`；已在 PR !150 README 中，属高级关键词 |
| `xc_parser` / `parse_xc` | PR !150 | 低 | `dft.md` 仅简要提及，待 !150 合并后完善 |
| `use_dm_only` 关键词 | PR !150 | 低 | 已在 PR !150 README 中，属高级关键词 |
| RRS-PBC 用户文档 | PR !150 | 低 | 已在 PR !150 README 中，专用功能 |
| 英文文档同步 | — | 中 | `solvent.md` 和 `geomopt.md` 需创建英文对应版本；`install.md` 英文版需同步更新 |
| Renormalized Doubles (FEAST) | #132 | 低 | 开发者文档候选 |
| RSH 解析梯度 | #111 | 低 | RSH 泛函梯度已自动生效，无需额外用户设置；可在 `dft.md` 中提及支持的梯度方法 |
| chkfile 基组投影 | PR !150 | 低 | `guessfile`/`chkfile` 基组不一致时自动投影的行为，已在 PR !150 README 中说明 |
