# Smelt 项目评估报告

## 结论

`Smelt` 当前已经不是“设计中”项目，而是一个已完成首版实现并可运行的 WGBS 下游分析 CLI 工具。代码结构清晰、模块边界明确、测试可通过，具备继续迭代的基础。

当前最主要的问题不是“没有实现”，而是对外文档不足和设计文档与实际实现存在若干偏差。这些偏差不会阻止项目运行，但会影响结果解释、用户预期和后续扩展的准确性。

## 项目概况

- 项目定位：WGBS 甲基化下游分析 CLI，覆盖 `site`、`window`、`element`、`metaplot`、`custom`、`dmr`、`stats` 七个子命令。
- 技术栈：Python + Typer + pandas/numpy/scipy + pyfaidx/pysam。
- 入口与组织方式：CLI 入口在 [smelt/cli.py](/home/hermes/Smelt/smelt/cli.py:1)，核心模块位于 [smelt/](/home/hermes/Smelt/smelt)。
- 设计说明主要在 [CLAUDE.md](/home/hermes/Smelt/CLAUDE.md:1)、[设计文档](/home/hermes/Smelt/docs/superpowers/specs/2026-04-29-smelt-design.md:1)、[实施计划](/home/hermes/Smelt/docs/superpowers/plans/2026-04-29-smelt-implementation.md:1)。

## 当前状态评估

- 代码实现完整度较高，七个主功能模块均已存在。
- 测试状态良好：本地执行 `pytest -q`，结果为 `81 passed, 14 warnings`。
- 仓库包含单元测试、CLI 测试和集成测试数据，说明作者已经考虑了基本验证路径。
- [README.md](/home/hermes/Smelt/README.md:1) 基本为空，导致项目对外可用性明显弱于内部实现成熟度。

## 优点

- 模块拆分合理，基本符合设计文档中的职责划分。
- 坐标转换规则明确，内部统一使用 0-based half-open，相关工具在 [smelt/utils.py](/home/hermes/Smelt/smelt/utils.py:1)。
- CLI 接口完整，参数设计与设计稿总体一致，见 [smelt/cli.py](/home/hermes/Smelt/smelt/cli.py:1)。
- I/O、统计、窗口聚合、DMR 调用等核心流程都已有实现，见 [smelt/io.py](/home/hermes/Smelt/smelt/io.py:1)、[smelt/window.py](/home/hermes/Smelt/smelt/window.py:1)、[smelt/dmr.py](/home/hermes/Smelt/smelt/dmr.py:1)。
- 测试覆盖范围比一般原型项目更好，说明代码不是一次性脚本式实现。

## 主要问题与风险

### 1. 文档层面

- [README.md](/home/hermes/Smelt/README.md:1) 未承担项目入口文档职责。
- 当前用户若不阅读内部文档，很难快速知道安装方式、输入格式、典型命令和输出示例。
- 这会直接影响项目可交付性和可复用性。

### 2. 设计与实现不完全一致

- 设计要求“严格染色体命名一致，发现不匹配时报错”，但当前实现中未看到统一的强校验逻辑，尤其在 `element`、`custom`、`metaplot` 路径上更明显，见 [smelt/element.py](/home/hermes/Smelt/smelt/element.py:1)、[smelt/custom.py](/home/hermes/Smelt/smelt/custom.py:1)、[smelt/metaplot.py](/home/hermes/Smelt/smelt/metaplot.py:1)。
- 设计要求 `element` 合并同一 feature 的连续片段，但当前实现是逐行 GTF 记录单独统计，没有做 feature-level merge，见 [smelt/element.py](/home/hermes/Smelt/smelt/element.py:1)。
- 设计要求 `metaplot` 输出包含 `sample` 维度，但当前实现输出只有 `region/bin/ratio` 聚合结果，见 [smelt/metaplot.py](/home/hermes/Smelt/smelt/metaplot.py:1)。
- 设计要求 `stats` 提供 genome-wide 和 per-chromosome 统计，当前实现只看到按 `chr, context` 聚合，没有 genome-wide 汇总行，见 [smelt/stats.py](/home/hermes/Smelt/smelt/stats.py:1)。

### 3. 生物学语义可能存在偏差

- `site` 计算中默认把所有位点 `strand` 设为 `"+"`，见 [smelt/site.py](/home/hermes/Smelt/smelt/site.py:1)。这与设计里对 CHH/CHG 链特异性的要求不完全一致。
- `FASTA` 上下文分类逻辑主要依赖当前位置下游序列，见 [smelt/site.py](/home/hermes/Smelt/smelt/site.py:1) 和 [smelt/utils.py](/home/hermes/Smelt/smelt/utils.py:1)。如果没有可靠链信息，这种实现更接近简化版上下文判定。
- `merge_cpg_strands()` 目前按同一 `chr,pos` 聚合，见 [smelt/utils.py](/home/hermes/Smelt/smelt/utils.py:1)。这未必等价于真正的 CpG 双链 dyad 合并，因为互补链 CpG 常常是相邻坐标而非同一坐标。

### 4. DMR 逻辑实现较简化

- 设计文档提到一些更细的排除场景，如“一组总覆盖为 0、另一组大于 0 时标记 low-coverage”，当前实现没有完整体现这类分类，见 [smelt/dmr.py](/home/hermes/Smelt/smelt/dmr.py:1)。
- 当前 `excluded` 输出里大量“未过阈值”的窗口也会被写入排除表，这在实用上可行，但语义上与“因数据问题被排除”混在一起。

## 测试评估

- 测试总体是项目的强项，覆盖了 CLI、I/O、窗口、DMR、统计与集成流程。
- 但测试更偏向“功能可运行”和“接口存在”，对上面这些设计偏差的约束还不够强。
- 换句话说，测试能证明当前实现稳定，但不能完全证明它与设计目标完全一致。

## 综合判断

`Smelt` 当前适合定义为：

“已完成 MVP/首版实现、测试通过、具备继续增强条件的 WGBS 分析工具”

它的主要短板在于：

- 对外文档缺失
- 若干设计承诺尚未严格落地
- 某些生物学细节实现偏简化

## 建议优先级

1. 先补齐 [README.md](/home/hermes/Smelt/README.md:1)，把安装、输入格式、命令示例、输出说明补完整。
2. 明确修正设计与实现的差异，优先处理 `element` feature merge、`stats` genome-wide 汇总、染色体命名校验。
3. 重新审视 `site` 的链信息与 CpG 双链合并逻辑，避免结果语义和设计稿不一致。
4. 为上述关键行为补充更严格的测试，而不只是保证“命令能跑通”。
