# Logger 日志系统

> **实验性功能**。Logger 日志系统在 [PR#175](https://gitee.com/restgroup/rest/pulls/175) 中引入，基于 [`log`](https://crates.io/crates/log) 和 [`env_logger`](https://crates.io/crates/env_logger) crate 构建。该系统的 API 和行为在未来版本中可能发生变化。

## 用法

### 输入卡关键词

`print_level` 关键词控制日志输出等级，其映射关系由 `utilities/log.rs` 中的 `printlevel2loglevel()` 函数定义：

| `print_level` | 日志等级 | 说明 |
|---|---|---|
| 0 | `Info` | 仅输出 info、warn、error 级别信息 |
| 1 | `Info` | 缺省值，同 0 |
| 2 | `Debug` | 额外输出 debug 级别调试信息 |
| ≥3 | `Trace` | 输出所有追踪信息 |

## 输出格式

日志系统使用自定义格式化输出：

- **`Info` 级别**：仅输出消息正文，无前缀
- **其他级别**（`Debug`、`Warn`、`Trace`、`Error`）：输出格式为 `[LEVEL module] message`

其中 `module` 为 Rust 模块路径的最后一段（例如 `scf_io::module` 显示为 `module`）。

## 实现

### 初始化

日志系统在 `ctrl_io/mod.rs` 的 `parse_ctrl_keywords()` 函数中初始化：

```rust
let mut builder = Builder::new();
builder.target(Target::Stdout);
builder.filter_level(printlevel2loglevel(tmp_input.print_level));
builder.format(|buf, record| {
    if record.level() == Level::Info {
        writeln!(buf, "{}", record.args())
    } else {
        writeln!(buf, "[{:<5} {}] {}", record.level(),
            record.target().split_once("::").map(|(_, rest)| rest).unwrap_or(record.target()),
            record.args())
    }
});
let _ = builder.try_init();
```

关键点：
- 日志输出目标为 stdout
- 使用 `try_init()` 而非 `init()`，以便在测试等场景中重复初始化时不会 panic

### 核心函数

`utilities/log.rs` 提供 `print_level` 到 `log::LevelFilter` 的映射：

```rust
pub fn printlevel2loglevel(level: usize) -> LevelFilter {
    match level {
        0 => LevelFilter::Info,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        3.. => LevelFilter::Trace,
        _ => LevelFilter::Info,
    }
}
```

### 运行时日志等级调整

在某些计算密集区域，可临时降低日志等级以减少输出开销。例如 `initial_guess/mod.rs` 在 SAD 初始猜测计算期间将日志等级临时设为 `Info`：

```rust
let cur_log_level = log::max_level();
log::set_max_level(LevelFilter::Info);
scf_data.density_matrix = initial_guess_from_sad(&scf_data.mol, mpi_operator);
log::set_max_level(cur_log_level);
```

## 使用日志宏的模块

- `ctrl_io/mod.rs` — `info!`、`debug!`、`warn!`
- `ctrl_io/path_util.rs` — `warn!`、`debug!`
- `scf_io/mod.rs` — `info!`、`debug!`、`trace!`、`warn!`
- `scf_io/scfrecord.rs` — `debug!`
- `ri_jk/pure_direct.rs` — `warn!`
- `ri_jk/pure_incore.rs` — `warn!`
- `initial_guess/mod.rs` — 运行时日志等级调整

## 注意事项

- 该功能为实验性，未来可能发生 API 或行为变更
- `log` 和 `env_logger` 为必需依赖，始终编译，无 feature flag 控制
- 建议新的模块使用 `info!`/`debug!`/`warn!`/`trace!` 宏代替 `println!`，以实现统一的输出控制
