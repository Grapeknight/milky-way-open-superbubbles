<a id="top"></a>

# The Milky Way as a "Phantom Galaxy": A Bubble-Dominated Disk Sculpted by Stellar Feedback and the Origin of the Radcliffe Wave

**[中文说明](#chinese) · [English documentation](#english)**

**代码 / Code:** [GitHub](https://github.com/Grapeknight/milky-way-open-superbubbles) · **完整代码、数据与结果 / Complete code, data and results:** [Releases v1.0.1](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/tag/v1.0.1)

GitHub 根目录直接存放完整包 `code/` 中的内容；完整 `code/data/results` 请从 Releases 下载。
The GitHub root contains the contents of the complete package's `code/` directory; download the complete `code/data/results` package from Releases.

<a id="chinese"></a>

## 中文说明

### 银河系如“幽灵星系”：恒星反馈塑造的气泡主导盘与拉德克利夫波的起源

本仓库提供论文 *The Milky Way as a "Phantom Galaxy": A Bubble-Dominated Disk Sculpted by Stellar Feedback and the Origin of the Radcliffe Wave* 的分析代码。配合 Releases 中的完整包，可复现研究中的数据分析、验证分析、图件、数据表和交互式可视化。

AI 辅助使用声明：整理、撰写说明和一致性检查过程中使用了 AI 工具。科学分析、数据解释和最终发布内容由作者负责。

与本 README 同目录的静态图 `pipeline_runtime_summary.png` 汇总了各编号脚本的参考实际运行时间，其中脚本 0b 的耗时是外部原始数据导出估计值。

完整包遵循 Code Ocean 标准布局：脚本和文档位于 `code/`，输入位于同级的 `data/`，输出位于同级的 `results/`。从 `code/` 工作目录执行时，输入路径为 `../data/`，中间产物路径为 `../results/intermediate_output/`，最终输出路径为 `../results/`。

### 下载完整发表包

[GitHub 代码仓库](https://github.com/Grapeknight/milky-way-open-superbubbles)的根目录直接放置完整包 `code/` 中的文件，包括 Python 脚本、本 README、文档和运行入口；仓库内没有额外的 `code/` 子目录。直接克隆仓库仅获得代码，不会下载 `data/` 和 `results/`。完整输入数据和已保存结果通过 [Releases v1.0.1](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/tag/v1.0.1) 分发，v1.0.0 保留为历史版本。GitHub 自动生成的 **Source code** 压缩包也只包含代码。

请下载以下三个 ZIP 文件：

- [milky-way-open-superbubbles-v1.0.1-code-data.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-code-data.zip)
- [milky-way-open-superbubbles-v1.0.1-results-1.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-results-1.zip)
- [milky-way-open-superbubbles-v1.0.1-results-2.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-results-2.zip)

将三个压缩包解压到同一个父目录，合并共同的 `milky-way-open-superbubbles-v1.0.1/` 文件夹，即可得到 `code/`、`data/`、`results/` 三个目录。完整包的说明仍位于 `code/README.md`。这些是互相独立的普通 ZIP 文件，不需要分卷解压工具。

运行时进入完整包的 `code/`。如使用 Git 克隆代码，运行环境仍需保持代码目录与 `data/`、`results/` 同级；GitHub 根目录的文件对应完整包的 `code/` 内容。不要向已经解压、非空的 `code/` 目录直接执行 `git clone`。

`SHA256SUMS.txt` 用于校验下载附件；`file-manifest.csv` 记录完整包内每个文件的大小及 SHA-256 校验值。Linux/macOS/WSL 可运行 `sha256sum -c SHA256SUMS.txt`（macOS 也可用 `shasum -a 256 -c SHA256SUMS.txt`）；PowerShell 可用 `Get-FileHash -Algorithm SHA256 <filename>`，与对应校验值比较。

### 输入数据及来源

主要输入星表和数据产品如下：

- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 15：全天三维尘埃图，通过 `dustmaps3d` 访问。
- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 16：用于追踪壳层的分子云星表。
- Wang et al. 2026, submitted：独立伴随工作提供的银河系尘埃泡星表 `../data/Bubbles.csv`。直接读取此表的是脚本 4、10、12、16、22、25；更新星表后应重新运行这六个脚本，使图件、表格和标注保持一致。
- Reid et al. 2019, ApJ, 885, 131：旋臂轨迹与大质量恒星形成区（HMSFR）样本。
- Hunt, E. L., and Reffert, S. 2024, Astron. Astrophys., 686, A42：疏散星团星表；相对脚本工作目录的路径为 `../data/star_cluster_data/hunt2024_clusters_full.csv`。
- Konietzka, R., et al. 2024, Nature, 628, 62–65：拉德克利夫波年轻星团样本 `../data/star_cluster_data/Konietzka2023.csv`。
- MWISP / Galactic Plane Survey `12CO`、CfA `12CO` 和 HI4PI HI 数据产品：紧凑 PV 缓存位于 `../data/COdata/` 和 `../data/HIdata/`。

### 脚本指南与运行时间

全部 27 个编号脚本构成完整分析流程。脚本 0b 用于原始 CO/HI 数据导出，需要外部原始巡天立方体；完整包已经包含其生成的紧凑 PV 产品，可直接供后续步骤复现使用。下表提供简要说明，每个脚本开头都有更详细的注释。

表中路径均相对 `code/`：输入为 `../data/...`，中间产物为 `../results/intermediate_output/...`，最终结果为 `../results/...`。耗时为记录的代表性实际运行时间，会随环境而变化；脚本 0b 为外部原始数据导出估计时间，因为原始巡天数据未包含在包中。

| 顺序 | 脚本 | 参考耗时 | 主要任务 | 主要输出及论文用途 |
|---:|---|---:|---|---|
| 0a | `0_Data_preparation_1.py` | 20 min 58.49 s | 使用 dustmaps3d，以限制内存占用的分块方式建立原始/平滑三维尘埃立方体、XY 平均尘埃图和拉德克利夫波尘埃截面。 | 在 `../results/intermediate_output/3d_dust_map_products/` 生成可复现中间数据，供脚本 5、8、11、12、23、24、25 使用。 |
| 0b | `0_Data_preparation_2.py` | > 6 h (外部原始数据导出) | 从外部原始 FITS 立方体导出紧凑的 MWISP/CfA/HI4PI 位置—速度（PV）缓存。数百 GB 的原始数据未随包分发，仅靠完整包不能运行此步骤。 | `../data/COdata/`、`../data/HIdata/`；已包含这些紧凑缓存，供脚本 21、22 使用。 |
| 1 | `1_OSBs_first_contact_MCs.py` | 3.27 s | 从每个开放超泡的初始位置沿射线识别首次接触的分子云。 | `../results/intermediate_output/1_first_contact_molecular_clouds/`；供脚本 2、3 使用。 |
| 2 | `2_OSBs_3D_fitting.py` | 0.91 s | 根据首次接触分子云拟合椭球帽或椭圆柱壳层几何。 | `../results/intermediate_output/2_automated_fit_results/`；供脚本 3、4 使用。 |
| 3 | `3_OSBs_MCMC_and_parameter_table.py` | 1 min 16.63 s | 利用 MCMC 估计不确定度，汇总 30 个目标的最终参数表。 | `../results/superbubble_final_fit_parameters.csv`、`../results/figures/3_mcmc_analysis.png` 和 corner 图；提供完整星表和 MCMC 诊断。 |
| 4 | `4_OSBs_plotting.py` | 1 h 02 min 35.22 s | 为每个超泡绘制六面板识别图。 | `../results/SB_figures/SB{N}.png`；用于补充信息和扩展数据中的超泡图。 |
| 5 | `5_OSBs_parameter_statistics.py` | 9 min 15.09 s | 测量空腔、壳层、密度脊统计量及证据分级参数。 | `../results/figures/5_parameter_histogram_matrix.png`；扩展数据参数分布图。 |
| 6 | `6_OSBs_initial_sensitivity.py` | 37 min 30.93 s | 扰动初始几何，重新检验首次接触分子云的恢复情况。 | `../results/figures/6_initial_sensitivity_overlap.png`；初始条件敏感性图。 |
| 7 | `7_OSBs_recall_analysis.py` | 0.95 s | 分析独立重复识别结果及不同识别者的一致性。 | `../results/figures/7_independent_repeat_recall.png`；独立识别召回率图。 |
| 8 | `8_OSBs_random_dust_test.py` | 3 h 12 min 21.40 s | 在固定拟合几何下运行三维 IAAFT 随机密度涨落检验。 | `../results/figures/8_condition_test_joint_p_all30.png`；随机密度检验图。 |
| 9 | `9_OSBs_shell_amplitude_test.py` | 35.22 s | 比较壳层密度跃变幅度与消光图、恒星样本的不确定度；需要时自动定位已安装的 dustmaps3d 缓存。 | `../results/figures/9_shell_amplitude_vs_extinction_error_3sigma.png`；壳层显著性图。 |
| 10 | `10_OSBs_bubbles_HMSFR_test.py` | 11 min 09.93 s | 统计与壳层关联的 A/B 级尘埃泡及大质量恒星形成区（HMSFR），并进行随机几何对照。 | `../results/figures/10_shell_assignment_xy.png`、`10_shell_count_random_test_overview.png`；恒星反馈示踪体关联图。 |
| 11 | `11_OSBs_top_view.py` | 21.01 s | 绘制银河系俯视开放超泡概览及 RGB 分层尘埃图。 | `../results/figures/11_topview_rgb_dust_layers.png`；正文图 1 的源面板。 |
| 12 | `12_OSBs_dust_slices.py` | 11.47 s | 绘制三个垂直尘埃截面，叠加开放超泡椭圆、尘埃泡和 HMSFR。 | `../results/figures/12_three_dust_slices.png`；正文图 2。 |
| 13 | `13_YSO_traceback_inputs.py` | 3.65 s | 建立 G1/G2 的 Hunt 与 Konietzka 年轻星团轨道回溯输入表。 | `../results/intermediate_output/13_traceback_cluster_input_data/`；供脚本 14、15、16、23、24 使用。 |
| 14 | `14_G1_YSO_traceback.py` | 3 min 24.56 s | 积分 G1 的 55 Myr 轨道，采用年龄加权的 HDBSCAN 恢复星团家族。 | `../results/intermediate_output/14_g1_traceback_clustering/`；供脚本 16、17 使用。 |
| 15 | `15_G2_YSO_traceback.py` | 1 min 26.16 s | 积分 G2 的 50 Myr 轨道，采用不加年龄权重的 HDBSCAN 恢复星团家族。 | `../results/intermediate_output/15_g2_traceback_clustering/`；供脚本 16 使用。 |
| 16 | `16_G1G2_YSO_traceback_plots.py` | 11.83 s | 绘制 G1/G2 回溯截面、尺度—年龄关系、残差速度场及家族—超泡组合面板，包含从 ../data/Bubbles.csv 读取的 HBW35 叠加轮廓。 | `../results/figures/traceback_cluster_figures/`；扩展数据回溯面板及年轻星团图。 |
| 17 | `17_SB31_RW_boundary_YSO_traceback.py` | 1.96 s | 为 G1 的 02 组及拉德克利夫波年轻星团绘制 SB31/拉德克利夫波边界回溯补充图。 | `../results/figures/traceback_cluster_figures/17_g1_family2_median_subtracted_traceback_slices.png`；SB31/拉德克利夫波回溯图。 |
| 18 | `18_YSO_OSBs_shell_random_test.py` | 25.27 s | 通过银经随机化检验不同年龄窗口的年轻星团在开放超泡内部和尘埃壳层上的分布。 | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png`；扩展数据年龄分层图。 |
| 19 | `19_OSBs_energy_budget.py` | 10.28 s | 采用蒙特卡洛能量预算模型估计动力学年龄及所需超新星数。 | `../results/figures/19_superbubble_energy_budget_age_sn.png`；能量预算图，并为脚本 20 提供输入。 |
| 20 | `20_OSBs_YSO_SN_estimation.py` | 1 min 37.11 s (首次运行，包含缓存构建) | 依据成员星测光约束，为每个星团匹配 100 组随机 IMF 星族；结合 PARSEC 等时线年龄及完整六维相空间轨道，筛选穿过各超泡随时间演化、按收缩关系校准体积的前身星，估计历史超新星供给；与面密度 922 和 2000–4000 Msun Myr^-1 kpc^-2 对应的局地恒星形成率供给比较。 | `../results/Open_superbubbles.csv`、`../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`；成员星贡献、各超泡超新星供给、样本筛选统计及方法参数。 |
| 21 | `21_OSBs_gas_PV_diagrams.py` | 27.98 s | 针对三种壳层速度模型绘制 CO/HI 银经—速度对比。 | `../results/figures/21_PV_*.png`；CO/HI PV 图。 |
| 22 | `22_bubbles_gas_PV_diagram.py` | 5.61 s | 在 MWISP CO PV 图上叠加低银纬 A/B 级尘埃泡。 | `../results/figures/22_gradeAB_low_lat_bubbles_PV_overlay.png`；闭合尘埃泡 PV 对比图。 |
| 23 | `23_RW_cross_sections.py` | 17.46 s | 绘制拉德克利夫波尘埃和年轻星团截面，叠加开放超泡壳层。 | `../results/figures/23_fig3.png`；正文图 3 的尘埃与年轻星团组合图。 |
| 24 | `24_RW_like_wave.py` | 12.18 s | 绘制另外三个类似拉德克利夫波的尘埃/星团截面。 | `../results/figures/24_rw_like_three_slices.png`；正文图 4。 |
| 25 | `25_interactive_three_volume_viewer.py` | 7.05 s | 生成嵌入科学数据的 Three.js 三维尘埃体可视化页面。 | `../results/interactive_3d_viewer.html`；开放超泡三维可视化的本地源文件。 |

### 论文图号与源文件对应

#### 正文图

| 正文图号 | 完整包中的源文件 | 对应脚本 |
|---|---|---|
| 正文图 1 | 组合图：银河系开放超泡/RGB 面板来自 `../results/figures/11_topview_rgb_dust_layers.png`；JWST/NGC 628 面板来自外部图像。 | 脚本 11 加外部 JWST 图像拼版 |
| 正文图 2 | `../results/figures/12_three_dust_slices.png` | 脚本 12 |
| 正文图 3 | `../results/figures/23_fig3.png` | 脚本 23 |
| 正文图 4 | `../results/figures/24_rw_like_three_slices.png` | 脚本 24 |

#### 扩展数据图

当前论文 PDF 中，扩展数据图延续四张正文图的编号，因此标为 Extended Data Fig. 5 至 Extended Data Fig. 12。

| 扩展数据图号 | 完整包中的源文件 | 对应脚本 |
|---|---|---|
| 扩展数据图 5 | `../results/figures/5_parameter_histogram_matrix.png` | 脚本 5 |
| 扩展数据图 6 | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png` | 脚本 18 |
| 扩展数据图 7 | `../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`；使用脚本 19 生成的能量预算汇总。 | 脚本 20，依赖脚本 19 的输出 |
| 扩展数据图 8 | `../results/figures/traceback_cluster_figures/G1G2_SBs_YSO.png` | 脚本 16 |
| 扩展数据图 9 | `../results/SB_figures/SB2.png` | 脚本 4 |
| 扩展数据图 10 | `../results/SB_figures/SB15.png`、`../results/SB_figures/SB16.png` | 脚本 4 |
| 扩展数据图 11 | `../results/SB_figures/SB25.png` | 脚本 4 |
| 扩展数据图 12 | `../results/SB_figures/SB31.png` | 脚本 4 |

扩展数据星表以脚本 3 生成的 `../results/superbubble_final_fit_parameters.csv` 和脚本 20 生成的简化星表 `../results/Open_superbubbles.csv` 为基础。

#### 补充信息

此处不逐项列出补充信息图号。补充信息中的所有开放超泡识别图均由 `4_OSBs_plotting.py` 生成，并保存在 `../results/SB_figures/`。

### 代码仓库与完整包结构

GitHub 代码仓库的根目录结构如下：

```text
milky-way-open-superbubbles/
|-- 0_Data_preparation_1.py
|-- 0_Data_preparation_2.py
|-- 1_OSBs_first_contact_MCs.py ... 25_interactive_three_volume_viewer.py
|-- final_data_products_description.md
|-- pipeline_runtime_summary.png
|-- README.md
|-- run
`-- run.sh
```

Releases 完整包保留 Code Ocean 的三个目录结构：

```text
milky-way-open-superbubbles-v1.0.1/
|-- code/
|   |-- 0_Data_preparation_1.py
|   |-- 0_Data_preparation_2.py
|   |-- 1_OSBs_first_contact_MCs.py ... 25_interactive_three_volume_viewer.py
|   |-- final_data_products_description.md
|   |-- pipeline_runtime_summary.png
|   |-- README.md
|   |-- run
|   `-- run.sh
|-- data/
|   |-- COdata/
|   |-- HIdata/
|   |-- Cantilever/
|   |-- star_cluster_data/
|   |-- Bubbles.csv
|   |-- MCs.csv
|   |-- OSB_initial_pre_identification_parameters.csv
|   |-- publication_plotting_parameters.csv
|   |-- Reid2019_spirals_xy.csv
|   |-- extinction_error.fits
|   `-- third_party_recall_test.csv
`-- results/
    |-- intermediate_output/       # 可复现的中间数据产品
    |-- figures/                   # 生成的论文图和诊断图
    |-- SB_figures/                # 生成的 OSB 六面板图
    |-- Open_superbubbles.csv
    |-- superbubble_final_fit_parameters.csv
    `-- interactive_3d_viewer.html
```


当前完整包的近似大小如下：

| 目录 | 近似大小 | 用途 |
|---|---:|---|
| `code/` | 1.5 MB | 脚本、说明、运行时间汇总图和 Code Ocean 运行入口 |
| `data/` | 0.64 GB | 输入星表和紧凑巡天数据产品 |
| `results/` | 2.58 GB | 最终结果和科学中间产物；脚本 20 运行时会重建其 IMF/轨道缓存 |

### 系统要求

- 操作系统：Windows 11、macOS 或 Linux。
- Python：建议 3.11，以匹配记录的参考环境。
- 磁盘：随包输入约 0.64 GB；完整流程及生成缓存建议至少预留 10 GB 可用空间。
- 网络：安装依赖、首次下载所需的 `dustmaps3d` 缓存，以及加载交互式查看器的 Three.js 模块时需要网络。
- 内存：完整流程建议至少 16 GB；`0_Data_preparation_2.py` 预处理大型原始 FITS 时，可能需要更多内存，具体取决于输入立方体。

Windows 说明：`4_OSBs_plotting.py` 的全天图需要 `healpy`。在本文记录的原生 Windows Python 环境下，请通过 WSL 运行流程。Code Ocean、Linux 或 macOS 安装 `healpy` 后可以直接运行。

### 已检查的参考环境

本 README 记录的当前本地环境如下：

| 组件 | 版本或说明 |
|---|---|
| Python | 3.11.3 |
| NumPy | 1.24.3 |
| pandas | 2.3.1 |
| Matplotlib | 3.7.2 |
| SciPy | 1.10.1 |
| Astropy | 7.1.0 |
| dustmaps3d | 2.2.13 |
| PyArrow | 20.0.0 |
| emcee | 3.1.4 |
| corner | 2.2.2 |
| galpy | 1.11.2 |
| hdbscan | 0.8.40 |
| astropy-healpix | 1.0.3 |
| Pillow | 9.4.0 |
| healpy | 当前原生 Windows Python 未安装；脚本 4 请使用 WSL/Linux |

以上是参考版本；其他环境可能产生细微数值或绘图差异。下面的安装命令直接列出所需 Python 包。

### 安装

解压完整包后，在包含 `code/`、`data/`、`results/` 的完整包根目录创建并激活虚拟环境。Linux、macOS 或 WSL 的完整流程环境可按以下命令安装：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas scipy matplotlib emcee corner pyarrow astropy galpy numba hdbscan healpy astropy-healpix dustmaps3d platformdirs Pillow tqdm
```

依赖中的 `numba` 供脚本 20 使用，`tqdm` 供可选的原始巡天导出进度显示使用。Windows 用户应在 WSL 内创建环境，使 `run.sh` 使用 Linux Python 和其中安装的 `healpy`。

首次使用 `dustmaps3d` 时可能下载 Wang 等人的三维尘埃图缓存。多数脚本所需的处理后尘埃数据，由 `0_Data_preparation_1.py` 生成并写入运行输出目录 `../results/intermediate_output/3d_dust_map_products/`。

### 快速开始

Code Ocean 中，将 `/code/run` 设为运行文件。入口会建立输出目录，并依次执行 26 个可在完整包基础上复现的脚本。

在本地 Linux 或 WSL 中，从完整包根目录运行：

```bash
cd code
bash run.sh
```

Windows 用户请进入 WSL，切换至完整包根目录，激活上述 WSL 虚拟环境后执行相同命令。原生 Windows Python 虚拟环境与 WSL 环境相互独立。

例外是 `0_Data_preparation_2.py`：它需要原始 HI4PI、CfA CO 和 MWISP FITS 立方体，不能仅凭完整包运行。原始数据有数百 GB，不随包分发；包中包含由这些原始数据导出的紧凑 PV 缓存，以及导出这些缓存的代码。

若逐个手动运行脚本，先从完整包根目录建立输出目录：

```bash
cd code
mkdir -p ../results ../results/intermediate_output
```

然后按顺序执行：

```powershell
# 在 ../results/intermediate_output/3d_dust_map_products/ 生成可复现的三维尘埃数据产品。
python 0_Data_preparation_1.py

# 可选：仅用于从原始巡天数据重新导出。
# 需要外部原始 HI4PI、CfA CO 和 Galactic Plane Survey / MWISP 数据。
# 完整包的 ../data/COdata/ 和 ../data/HIdata/ 已包含后续所需的紧凑 PV 缓存。
# python 0_Data_preparation_2.py --mwisp-fits <path> --cfa-fits <path> --hi4pi-fits <path>

python 1_OSBs_first_contact_MCs.py
python 2_OSBs_3D_fitting.py
python 3_OSBs_MCMC_and_parameter_table.py
python 4_OSBs_plotting.py
python 5_OSBs_parameter_statistics.py
python 6_OSBs_initial_sensitivity.py
python 7_OSBs_recall_analysis.py
python 8_OSBs_random_dust_test.py
python 9_OSBs_shell_amplitude_test.py
python 10_OSBs_bubbles_HMSFR_test.py
python 11_OSBs_top_view.py
python 12_OSBs_dust_slices.py
python 13_YSO_traceback_inputs.py
python 14_G1_YSO_traceback.py
python 15_G2_YSO_traceback.py
python 16_G1G2_YSO_traceback_plots.py
python 17_SB31_RW_boundary_YSO_traceback.py
python 18_YSO_OSBs_shell_random_test.py
python 19_OSBs_energy_budget.py
python 20_OSBs_YSO_SN_estimation.py
python 21_OSBs_gas_PV_diagrams.py
python 22_bubbles_gas_PV_diagram.py
python 23_RW_cross_sections.py
python 24_RW_like_wave.py
python 25_interactive_three_volume_viewer.py
```

脚本 4 不能在本文记录的原生 Windows Python 环境中直接运行。Windows 请在 WSL 中执行 `bash run.sh`；Code Ocean、Linux、macOS 安装 `healpy` 后可直接执行 `python 4_OSBs_plotting.py`。

### 数据准备说明

两个 `0_` 脚本生成后续步骤共用的数据产品。

`0_Data_preparation_1.py` 生成：

- `../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/xy_mean_dust_map.parquet`
- `../results/intermediate_output/3d_dust_map_products/radcliffe_wave_dust_slice.csv`

`0_Data_preparation_2.py` 记录从原始巡天立方体导出紧凑 CO/HI PV 缓存的方法，其相对代码目录的输出路径为：

- `../data/COdata/MWISP_12CO_PV_cache.npz`
- `../data/COdata/CfA_12CO_PV_cache.npz`
- `../data/HIdata/HI4PI_HI_PV_cache.npy`
- `../data/HIdata/HI4PI_HI_PV_axes.npz`
- `../data/HIdata/HI4PI_HI_PV_cache_summary.json`

这些文件实际保存在完整包根目录下的 `data/COdata/` 和 `data/HIdata/` 中。

如果不另外提供原始 HI4PI、CfA CO 和 Galactic Plane Survey / MWISP FITS 立方体，就不能直接重跑此导出脚本。数百 GB 的原始数据未分发；脚本 21、22 所需的紧凑缓存已经包含。要从原始数据重建缓存，请显式传入路径：

```powershell
python 0_Data_preparation_2.py `
  --mwisp-fits <path-to-mosaic_12CO.fits> `
  --cfa-fits <path-to-CfA-FITS-cube> `
  --hi4pi-fits <path-to-HI4PI-CAR.fits>
```

MWISP FITS 路径也可通过环境变量 `MWISP_DATA_PATH` 指定。

### 本地交互式三维可视化

从完整包根目录执行以下命令生成本地页面：

```bash
cd code
python3 25_interactive_three_volume_viewer.py
```

主要输入：

- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../data/Bubbles.csv`
- `../results/superbubble_final_fit_parameters.csv`
- `../results/Open_superbubbles.csv`（可选的简化显示表；如果缺失，查看器会改用完整参数表）

输出：

- `../results/interactive_3d_viewer.html`

HTML 内嵌 uint8 尘埃体纹理、开放超泡几何和尘埃泡标注。它从 `unpkg.com` 加载 Three.js 0.160.0 及附加模块，因此查看时需要网络。可选参数 `--output-bin <path>` 会另行导出一份体数据载荷；HTML 本身不依赖该文件。在 Code Ocean 中，运行后可从 `/results/interactive_3d_viewer.html` 获取页面。

### 输出约定

- `../data/` 包含原始星表和可重复使用的非本流程生成缓存。
- `../results/intermediate_output/` 包含可复现产物，例如生成的三维尘埃数据、各脚本的 CSV/JSON 汇总、随机试验和诊断表。
- `../results/` 包含流程运行后的论文最终数据产品。
- `../results/superbubble_final_fit_parameters.csv` 为脚本 3 生成的 30 个目标的权威几何参数表。
- `../results/Open_superbubbles.csv` 为脚本 20 在全部上游数据就绪后生成的简化星表。几何和 MCMC 不确定度来自脚本 3，壳层和尘埃指标来自脚本 5，随机密度检验 p 值来自脚本 8，动力学年龄和所需超新星数来自脚本 19；发表的证据等级来自脚本 20 中编码的简化星表等级映射。
- 与本 README 同目录的 `final_data_products_description.md` 仅为说明文件，解释最终数据产品及各列含义，不是分析输出。
- 同目录的 `pipeline_runtime_summary.png` 为规划计算流程使用的静态耗时汇总图，不由流程重新生成。
- 脚本 4 运行后，`../results/SB_figures/` 包含全部 30 个目标的六面板图。
- `../results/figures/` 包含脚本 3 和 5–24 生成的正文图、扩展数据图及诊断图。

多数脚本会将工作目录设为其所在的 `code/` 目录，并使用相对代码目录的路径；Code Ocean 和本地完整包使用相同的直接相对路径。GitHub 仓库根目录即对应这里的 `code/`。

### 可复现性说明

脚本 8 在填补缺失值、按需进行块平均后，使用包含空腔、壳层及周围环境的完整局地原始三维尘埃立方体设定 IAAFT 的傅里叶振幅。它精确保留密度分布，并近似保留功率谱。默认迭代 8 次；迭代次数本身不能保证收敛。详见[脚本 8 内的方法及诊断说明](8_OSBs_random_dust_test.py)。

30 个目标分别在 `../results/figures/random_density_diagnostics/` 中保存 `8_SB{id}_original_vs_random.png`（仅 PNG），采用第一组随机实现和标准的 8 次迭代设置。每张对照图的前三行是匹配位置的 XY、XZ、YZ 截面，左列为原始数据，右列为随机实现；下方三个面板展示完整立方体的三维功率谱、功率比和密度分布。图中保留简洁标题、坐标轴、单位与图例；详细方法及诊断数值记入说明和元数据，供撰写发表图注时使用。

六张空间截面共享 `Spectral_r` 色标，在实际检验网格（可选重采样之后）上，分别对每个二维截面直接施加仅用于显示的高斯平滑，面内每个轴的 sigma 均为 1 个体素。功率谱、密度分布和保存数组均使用显示平滑前的立方体；细节见脚本内方法说明。数组、频谱表和元数据保存在 `../results/intermediate_output/8_random_density_fluctuation_test/diagnostics/`。诊断元数据中的路径相对 `code/`，便于迁移；路径约定不改变科学统计量。

在 `code/` 中运行 `python 8_OSBs_random_dust_test.py --targets all --diagnostics-only --iaaft-max-iter 8`，可为每个目标生成一张对照图，而不改写蒙特卡洛样本和概率表。独立诊断运行可通过 `--diagnostic-dir`、`--comparison-fig-dir` 指定输出位置。

- `0_Data_preparation_2.py` 处理原始 HI4PI、CfA CO 和 Galactic Plane Survey / MWISP FITS 数据；数百 GB 的原始输入未分发，紧凑的后续 PV 缓存已包含。
- 脚本 8 运行三维 IAAFT 随机场检验，是标准分析流程中最耗时的脚本。
- 脚本 9 自动检查是否安装 `dustmaps3d`，并在可用时自动定位本地三维尘埃图数据包。正常运行只需标准 `dustmaps3d` 缓存，无须额外指定数据路径。
- 脚本 14、15 默认可以复用轨道缓存；如需完整重算，使用 `--recompute-orbit --recompute-labels`，耗时可能增加。
- 正文图 1 按前文说明由生成的面板与外部图像组合；脚本 23 直接输出完整的正文图 3。
- 即使删除 `../results/intermediate_output/` 和 `../results/`，只要保留完整包的 `../data/` 且外部 `dustmaps3d` 缓存可用，按顺序运行脚本即可重建可复现产物。

### 配套在线三维尘埃图平台

开放超泡可视化的配套在线三维尘埃图平台：

- https://nadc.china-vo.org/data/dustmaps/OSBs

### 许可证

MIT。

### 联系方式

Haibo Yuan 教授（Prof. Haibo Yuan）：yuanhb@bnu.edu.cn

Tao Wang：wt@mail.bnu.edu.cn

[返回顶部 / Back to top](#top) · [English documentation](#english)

---

<a id="english"></a>

## English documentation

This repository provides the code for the end-to-end publication package for the paper *The Milky Way as
a "Phantom Galaxy": A Bubble-Dominated Disk Sculpted by Stellar Feedback and
the Origin of the Radcliffe Wave*. It reproduces the data analysis, validation
analyses, figures, data tables, and interactive visualization used in the study.

AI assistance disclosure: AI tools were used to help organize, document, and
consistency-check this code package. The scientific analysis, data
interpretation, and final package content remain the responsibility of the
authors.

The static runtime overview figure `pipeline_runtime_summary.png` in this
`code/` directory summarizes reference wall-clock runtimes of the numbered
scripts. Script 0b is an external raw-data export estimate.

This package is organized for Code Ocean's standard capsule layout. Scripts and
documentation are in `code/`, input data are in `data/` at the
capsule root (`../data/` from `code/`), and reproducible-run outputs are written
to `../results/` from the `code/` working directory.

Scripts run from `code/` and use direct relative paths in the current file tree:
input files are read from `../data/`, intermediate outputs are written to
`../results/intermediate_output/`, and final outputs are written to
`../results/`. 

### Download the complete publication package

The [GitHub repository](https://github.com/Grapeknight/milky-way-open-superbubbles)
contains the contents of the complete package's `code/` directory directly at
the repository root: Python scripts, this README, documentation, and run entry
points. There is no extra `code/` subdirectory in the GitHub repository.
Cloning the repository downloads code only; it does not download `data/` or
`results/`. Input data and saved results
are distributed as assets of [release v1.0.1](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/tag/v1.0.1).
GitHub's automatically generated **Source code** archives contain the code only.

Version v1.0.0 remains available as a historical release.
For the complete package, download these three ZIP files from release v1.0.1:

- [milky-way-open-superbubbles-v1.0.1-code-data.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-code-data.zip)
- [milky-way-open-superbubbles-v1.0.1-results-1.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-results-1.zip)
- [milky-way-open-superbubbles-v1.0.1-results-2.zip](https://github.com/Grapeknight/milky-way-open-superbubbles/releases/download/v1.0.1/milky-way-open-superbubbles-v1.0.1-results-2.zip)

Extract all three archives into the same parent directory, merging the common
`milky-way-open-superbubbles-v1.0.1/` folder. The resulting capsule contains
`code/`, `data/`, and `results/`; the main instructions remain in `code/README.md`.
The archives are ordinary independent ZIP files and do not require a split-archive utility.

Run scripts from the complete package's `code/` directory. If using a Git
checkout, keep the code directory alongside `data/` and `results/`; files at
the GitHub root correspond to the complete package's `code/` contents. Do not
run `git clone` into an already extracted, nonempty `code/` directory.

`SHA256SUMS.txt` verifies the downloaded assets, and `file-manifest.csv` records
the size and SHA-256 checksum of every file inside the complete capsule.
On Linux/macOS/WSL, verify downloads with `sha256sum -c SHA256SUMS.txt` (or
`shasum -a 256 -c SHA256SUMS.txt` on macOS). On PowerShell, compare
`Get-FileHash -Algorithm SHA256 <filename>` with the corresponding checksum.

### Input data and sources

The main input catalogues and data products are:

- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 15: all-sky 3D dust map,
  accessed through `dustmaps3d`.
- Wang, T., Yuan, H., et al. 2025, ApJS, 280, 16: molecular-cloud catalogue
  used for shell-tracing clouds.
- Wang et al. 2026, submitted: Galactic dust-bubble catalogue
  (`../data/Bubbles.csv`) from a separate companion work.
  The direct consumers are scripts 4, 10, 12, 16, 22, and 25; rerun all six
  after replacing this catalogue so their figures, tables, and annotations
  remain synchronized.
- Reid et al. 2019, ApJ, 885, 131: spiral-arm loci and HMSFR sample.
- Hunt, E. L., and Reffert, S. 2024, Astron. Astrophys., 686, A42:
  open-cluster catalogue in `../data/star_cluster_data/hunt2024_clusters_full.csv`
  from the script working directory.
- Konietzka, R., et al. 2024, Nature, 628, 62-65: Radcliffe-Wave young-cluster
  sample in `../data/star_cluster_data/Konietzka2023.csv`.
- MWISP / Galactic Plane Survey `12CO`, CfA `12CO`, and HI4PI HI products:
  compact PV caches under `../data/COdata/` and `../data/HIdata/`.

### Program guide and runtime

All numbered scripts are part of the complete analysis workflow. Script 0b is
the raw CO/HI export step and requires external raw survey cubes; the compact
PV outputs from that step are already included for downstream reproduction.
Shorter descriptions are given here; each script contains a more detailed header
at the top of the file.

In the table, paths are direct relative paths from `code/`: generated intermediate products are under `../results/intermediate_output/...`, generated final products are under `../results/...`, and inputs are under `../data/...`. Runtimes are representative recorded wall-clock measurements and depend on the execution environment. Script 0b is an external raw-data export estimate because the raw survey cubes are not packaged.

| Order | Script | Runtime | Main task | Main outputs and manuscript use |
|---:|---|---:|---|---|
| 0a | `0_Data_preparation_1.py` | 20 min 58.49 s | Build raw/smoothed 3D dust cubes, XY mean dust map, and Radcliffe-Wave dust slice from `dustmaps3d` using memory-bounded chunked output. | Reproducible intermediate inputs in `../results/intermediate_output/3d_dust_map_products/` for scripts 5, 8, 11, 12, 23, 24, 25. |
| 0b | `0_Data_preparation_2.py` | > 6 h (external raw-data export) | Export compact MWISP/CfA/HI4PI PV caches from external original FITS cubes. This step cannot be run from the package alone because the several-hundred-GB raw cubes are not included. | `../data/COdata/`, `../data/HIdata/`; these compact caches are included and used by scripts 21 and 22. |
| 1 | `1_OSBs_first_contact_MCs.py` | 3.27 s | Identify first-contact molecular clouds along rays for each open-superbubble seed. | `../results/intermediate_output/1_first_contact_molecular_clouds/`; input to scripts 2 and 3. |
| 2 | `2_OSBs_3D_fitting.py` | 0.91 s | Fit ellipsoidal-cap or elliptical-cylinder shell geometry to first-contact clouds. | `../results/intermediate_output/2_automated_fit_results/`; input to scripts 3 and 4. |
| 3 | `3_OSBs_MCMC_and_parameter_table.py` | 1 min 16.63 s | Estimate MCMC uncertainties and assemble the final 30-object parameter table. | `../results/superbubble_final_fit_parameters.csv`, `../results/figures/3_mcmc_analysis.png`, corner plots; source of the full catalogue table and MCMC diagnostics. |
| 4 | `4_OSBs_plotting.py` | 1 h 02 min 35.22 s | Draw the six-panel identification map for each superbubble. | `../results/SB_figures/SB{N}.png`; used in the Supplementary Information and Extended Data SB panels. |
| 5 | `5_OSBs_parameter_statistics.py` | 9 min 15.09 s | Measure cavity/shell/ridge statistics and evidence-grade parameters. | `../results/figures/5_parameter_histogram_matrix.png`; Extended Data parameter histogram. |
| 6 | `6_OSBs_initial_sensitivity.py` | 37 min 30.93 s | Perturb initial geometry and retest first-contact cloud recovery. | `../results/figures/6_initial_sensitivity_overlap.png`; initial-condition sensitivity figure. |
| 7 | `7_OSBs_recall_analysis.py` | 0.95 s | Analyze independent repeat-identification and observer agreement. | `../results/figures/7_independent_repeat_recall.png`; independent-identification recall figure. |
| 8 | `8_OSBs_random_dust_test.py` | 3 h 12 min 21.40 s | Run 3D IAAFT random-density-fluctuation tests at fixed fitted geometry. | `../results/figures/8_condition_test_joint_p_all30.png`; random-density test figure. |
| 9 | `9_OSBs_shell_amplitude_test.py` | 35.22 s | Compare shell jump amplitudes with extinction-map and stellar-sample uncertainties; automatically detects an installed `dustmaps3d` cache when needed. | `../results/figures/9_shell_amplitude_vs_extinction_error_3sigma.png`; shell-significance figure. |
| 10 | `10_OSBs_bubbles_HMSFR_test.py` | 11 min 09.93 s | Count shell-associated Grade A/B bubbles and HMSFRs; run random geometry controls. | `../results/figures/10_shell_assignment_xy.png`, `10_shell_count_random_test_overview.png`; feedback-tracer association figures. |
| 11 | `11_OSBs_top_view.py` | 21.01 s | Draw face-on open-superbubble overview and RGB layered dust map. | `../results/figures/11_topview_rgb_dust_layers.png`; source panel for main-text Fig. 1. |
| 12 | `12_OSBs_dust_slices.py` | 11.47 s | Draw three vertical dust slices with open-superbubble ellipses, bubbles, and HMSFRs. | `../results/figures/12_three_dust_slices.png`; main-text Fig. 2. |
| 13 | `13_YSO_traceback_inputs.py` | 3.65 s | Build G1/G2 Hunt and Konietzka young-cluster traceback input tables. | `../results/intermediate_output/13_traceback_cluster_input_data/`; input to scripts 14, 15, 16, 23, 24. |
| 14 | `14_G1_YSO_traceback.py` | 3 min 24.56 s | Integrate G1 orbits for 55 Myr and run age-weighted HDBSCAN family recovery. | `../results/intermediate_output/14_g1_traceback_clustering/`; input to script 16 and 17. |
| 15 | `15_G2_YSO_traceback.py` | 1 min 26.16 s | Integrate G2 orbits for 50 Myr and run no-age-weighted HDBSCAN family recovery. | `../results/intermediate_output/15_g2_traceback_clustering/`; input to script 16. |
| 16 | `16_G1G2_YSO_traceback_plots.py` | 11.83 s | Plot G1/G2 traceback slices, size-age curves, residual velocity fields, and the combined family-SB panel, including the HBW35 overlay read from `../data/Bubbles.csv`. | `../results/figures/traceback_cluster_figures/`; Extended Data traceback panel and young-cluster figures. |
| 17 | `17_SB31_RW_boundary_YSO_traceback.py` | 1.96 s | Produce the SB31/Radcliffe-Wave boundary traceback supplement for G1 group 02 and young RW clusters. | `../results/figures/traceback_cluster_figures/17_g1_family2_median_subtracted_traceback_slices.png`; SB31/RW traceback figure. |
| 18 | `18_YSO_OSBs_shell_random_test.py` | 25.27 s | Test young-cluster age windows inside open-superbubble volumes and on dust shells with longitude randomization. | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png`; Extended Data age-stratification figure. |
| 19 | `19_OSBs_energy_budget.py` | 10.28 s | Estimate dynamical ages and required supernova counts by Monte Carlo energy-budget modelling. | `../results/figures/19_superbubble_energy_budget_age_sn.png`; energy-budget figure and input to script 20. |
| 20 | `20_OSBs_YSO_SN_estimation.py` | 1 min 37.11 s (first run, including cache build) | Estimate historical cluster SNe by matching member-photometry constraints to 100 stochastic IMF populations per cluster, then gating sampled progenitors by PARSEC isochrone ages and complete-6D orbits through each bubble's evolving, contraction-calibrated volume; compare with local-SFR supply using 922 and 2000–4000 Msun Myr^-1 kpc^-2. | `../results/Open_superbubbles.csv`, `../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`; member-level contributions, per-bubble SN supply, selection statistics, and method parameters. |
| 21 | `21_OSBs_gas_PV_diagrams.py` | 27.98 s | Draw CO/HI longitude-velocity comparisons for three shell velocity prescriptions. | `../results/figures/21_PV_*.png`; CO/HI PV figures. |
| 22 | `22_bubbles_gas_PV_diagram.py` | 5.61 s | Overlay low-latitude Grade A/B bubbles on MWISP CO PV diagrams. | `../results/figures/22_gradeAB_low_lat_bubbles_PV_overlay.png`; closed-bubble PV comparison. |
| 23 | `23_RW_cross_sections.py` | 17.46 s | Draw Radcliffe-Wave dust and young-cluster cross sections with open-superbubble shell overlays. | `../results/figures/23_fig3.png`; combined dust and young-cluster panels for main-text Fig. 3. |
| 24 | `24_RW_like_wave.py` | 12.18 s | Draw three additional Radcliffe-Wave-like dust/cluster slices. | `../results/figures/24_rw_like_three_slices.png`; main-text Fig. 4. |
| 25 | `25_interactive_three_volume_viewer.py` | 7.05 s | Build the Three.js 3D dust-volume viewer with embedded scientific data. | `../results/interactive_3d_viewer.html`; local source for the OSB 3D visualization. |

### Figure provenance by document figure number

#### Manuscript figures

| Manuscript figure | Package source product | Corresponding script |
|---|---|---|
| Manuscript Fig. 1 | Composite figure. The Milky Way OSB/RGB panel comes from `../results/figures/11_topview_rgb_dust_layers.png`; the JWST/NGC 628 panel is external. | Script 11 plus external JWST-image assembly |
| Manuscript Fig. 2 | `../results/figures/12_three_dust_slices.png` | Script 12 |
| Manuscript Fig. 3 | `../results/figures/23_fig3.png` | Script 23 |
| Manuscript Fig. 4 | `../results/figures/24_rw_like_three_slices.png` | Script 24 |

#### Extended Data figures

In the current manuscript PDF, the Extended Data figures are printed as
Extended Data Fig. 5 through Extended Data Fig. 12 because the figure counter
continues after the four main-text figures.

| Extended Data figure | Package source product | Corresponding script |
|---|---|---|
| Extended Data Fig. 5 | `../results/figures/5_parameter_histogram_matrix.png` | Script 5 |
| Extended Data Fig. 6 | `../results/figures/18_hunt_young_clusters_shell_relative_offset_trend.png` | Script 18 |
| Extended Data Fig. 7 | `../results/figures/20_cluster_sn_estimate_in_superbubble_volume.png`; this figure uses the energy-budget summary generated by script 19. | Script 20, using script 19 outputs |
| Extended Data Fig. 8 | `../results/figures/traceback_cluster_figures/G1G2_SBs_YSO.png` | Script 16 |
| Extended Data Fig. 9 | `../results/SB_figures/SB2.png` | Script 4 |
| Extended Data Fig. 10 | `../results/SB_figures/SB15.png`, `../results/SB_figures/SB16.png` | Script 4 |
| Extended Data Fig. 11 | `../results/SB_figures/SB25.png` | Script 4 |
| Extended Data Fig. 12 | `../results/SB_figures/SB31.png` | Script 4 |

The Extended Data catalogue table is based on
`../results/superbubble_final_fit_parameters.csv` and the compact
`../results/Open_superbubbles.csv`, generated by scripts 3 and 20,
respectively.

#### Supplementary Information

The Supplementary Information figure provenance is not itemized here. All
Supplementary Information open-superbubble identification figures are generated
by `4_OSBs_plotting.py` and stored in `../results/SB_figures/`.

### Repository and complete-package contents

The GitHub code repository has this root layout:

```text
milky-way-open-superbubbles/
|-- 0_Data_preparation_1.py
|-- 0_Data_preparation_2.py
|-- 1_OSBs_first_contact_MCs.py ... 25_interactive_three_volume_viewer.py
|-- final_data_products_description.md
|-- pipeline_runtime_summary.png
|-- README.md
|-- run
`-- run.sh
```

The complete release package preserves the Code Ocean capsule layout:

```text
milky-way-open-superbubbles-v1.0.1/
|-- code/
|   |-- 0_Data_preparation_1.py
|   |-- 0_Data_preparation_2.py
|   |-- 1_OSBs_first_contact_MCs.py ... 25_interactive_three_volume_viewer.py
|   |-- final_data_products_description.md
|   |-- pipeline_runtime_summary.png
|   |-- README.md
|   |-- run
|   `-- run.sh
|-- data/
|   |-- COdata/
|   |-- HIdata/
|   |-- Cantilever/
|   |-- star_cluster_data/
|   |-- Bubbles.csv
|   |-- MCs.csv
|   |-- OSB_initial_pre_identification_parameters.csv
|   |-- publication_plotting_parameters.csv
|   |-- Reid2019_spirals_xy.csv
|   |-- extinction_error.fits
|   `-- third_party_recall_test.csv
`-- results/
    |-- intermediate_output/       # reproducible intermediate products
    |-- figures/                   # generated manuscript and diagnostic figures
    |-- SB_figures/                # generated six-panel OSB figures
    |-- Open_superbubbles.csv
    |-- superbubble_final_fit_parameters.csv
    `-- interactive_3d_viewer.html
```

Approximate package sizes in the current local copy:

| Directory | Approximate size | Role |
|---|---:|---|
| `code/` | 1.5 MB | Scripts, documentation, runtime summary, and Code Ocean run entry point |
| `data/` | 0.64 GB | Input catalogues and compact survey products |
| `results/` | 2.58 GB | Final outputs and scientific intermediate products; Script 20 rebuilds its IMF/orbit caches when run |

### System requirements

- Operating systems: Windows 11, macOS, or Linux.
- Python: 3.11 recommended to match the documented reference environment.
- Disk space: packaged inputs are about 0.64 GB; at least 10 GB free space is
  recommended for the full workflow and generated caches.
- Internet: required for dependency installation, the first `dustmaps3d` cache
  download if needed, and loading the interactive viewer's Three.js modules.
- Memory: at least 16 GB recommended for the full workflow; large raw FITS
  preprocessing in `0_Data_preparation_2.py` can require more, depending on the
  input cubes.

Windows note: `4_OSBs_plotting.py` uses `healpy` for all-sky panels. In the
tested native Windows Python environment, run the workflow through WSL. On
Code Ocean, Linux, or macOS, it can be run directly if `healpy` is installed.

### Tested environment

The current local environment checked for this README used:

| Component | Version or note |
|---|---|
| Python | 3.11.3 |
| NumPy | 1.24.3 |
| pandas | 2.3.1 |
| Matplotlib | 3.7.2 |
| SciPy | 1.10.1 |
| Astropy | 7.1.0 |
| dustmaps3d | 2.2.13 |
| PyArrow | 20.0.0 |
| emcee | 3.1.4 |
| corner | 2.2.2 |
| galpy | 1.11.2 |
| hdbscan | 0.8.40 |
| astropy-healpix | 1.0.3 |
| Pillow | 9.4.0 |
| healpy | Not installed in the native Windows Python; use WSL/Linux for script 4 |

The environment above records the reference versions; other environments can
produce small numerical or rendering differences. The install commands below
list the required Python packages directly.

### Installation

After extracting the complete release package, create and activate a virtual
environment from the capsule root (the directory containing `code/`, `data/`,
and `results/`). For the
complete workflow on Linux, macOS, or WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas scipy matplotlib emcee corner pyarrow astropy galpy numba hdbscan healpy astropy-healpix dustmaps3d platformdirs Pillow tqdm
```

The dependency list includes `numba` for Script 20 and `tqdm` for the optional
raw-survey export progress display. On Windows, create this environment inside
WSL so the Linux Python and `healpy` installation are used by `run.sh`.

First use of `dustmaps3d` may download the Wang et al. 3D dust-map cache. The
processed dust products needed by most scripts are generated by
`0_Data_preparation_1.py` under `../results/intermediate_output/3d_dust_map_products/`
inside the run output area.

### Quick start

On Code Ocean, set `/code/run` as the run file. The run script creates the
output directories and executes the 26 reproducible scripts in order.

For local Linux or WSL use from the capsule root:

```bash
cd code
bash run.sh
```

On Windows, open WSL, change to the complete capsule root, activate the WSL virtual
environment created above, and run the same commands. A native Windows Python
virtual environment is separate from the WSL environment.

The exception is `0_Data_preparation_2.py`: it cannot be executed from this
package alone because it requires the original HI4PI, CfA CO, and MWISP FITS
cubes. Those raw data are several hundred GB and are not included here. This
package includes the compact PV cache products and the code used to export those
PV products from the raw data.

For manual, individual-script execution, first prepare the output directories:

```bash
cd code
mkdir -p ../results ../results/intermediate_output
```

Then run the scripts in order:

```powershell
# Builds the reproducible 3D dust products under ../results/intermediate_output/3d_dust_map_products/.
python 0_Data_preparation_1.py

# Optional from-raw-data step only.
# Requires external raw HI4PI, CfA CO, and Galactic Plane Survey / MWISP data.
# The resulting compact PV caches are already included in ../data/COdata/ and ../data/HIdata/.
# python 0_Data_preparation_2.py --mwisp-fits <path> --cfa-fits <path> --hi4pi-fits <path>

python 1_OSBs_first_contact_MCs.py
python 2_OSBs_3D_fitting.py
python 3_OSBs_MCMC_and_parameter_table.py
python 4_OSBs_plotting.py
python 5_OSBs_parameter_statistics.py
python 6_OSBs_initial_sensitivity.py
python 7_OSBs_recall_analysis.py
python 8_OSBs_random_dust_test.py
python 9_OSBs_shell_amplitude_test.py
python 10_OSBs_bubbles_HMSFR_test.py
python 11_OSBs_top_view.py
python 12_OSBs_dust_slices.py
python 13_YSO_traceback_inputs.py
python 14_G1_YSO_traceback.py
python 15_G2_YSO_traceback.py
python 16_G1G2_YSO_traceback_plots.py
python 17_SB31_RW_boundary_YSO_traceback.py
python 18_YSO_OSBs_shell_random_test.py
python 19_OSBs_energy_budget.py
python 20_OSBs_YSO_SN_estimation.py
python 21_OSBs_gas_PV_diagrams.py
python 22_bubbles_gas_PV_diagram.py
python 23_RW_cross_sections.py
python 24_RW_like_wave.py
python 25_interactive_three_volume_viewer.py
```

Script 4 is not expected to run in the tested native Windows Python
environment. On Windows, run `bash run.sh` through WSL; on Code Ocean, Linux, or
macOS, run `python 4_OSBs_plotting.py` directly if `healpy` is installed.

### Data preparation notes

The two `0_` scripts generate shared data products used by later scripts.

`0_Data_preparation_1.py` builds:

- `../results/intermediate_output/3d_dust_map_products/raw_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../results/intermediate_output/3d_dust_map_products/xy_mean_dust_map.parquet`
- `../results/intermediate_output/3d_dust_map_products/radcliffe_wave_dust_slice.csv`

`0_Data_preparation_2.py` is provided to document how the compact CO/HI PV
caches were exported from the original raw survey cubes:

- code-relative `../data/COdata/MWISP_12CO_PV_cache.npz`
- code-relative `../data/COdata/CfA_12CO_PV_cache.npz`
- code-relative `../data/HIdata/HI4PI_HI_PV_cache.npy`
- code-relative `../data/HIdata/HI4PI_HI_PV_axes.npz`
- code-relative `../data/HIdata/HI4PI_HI_PV_cache_summary.json`

These files are physically stored under `../data/COdata/` and
`../data/HIdata/` in the capsule root.

This script cannot be run directly from the distributed package unless the user
also supplies the original HI4PI, CfA CO, and Galactic Plane Survey / MWISP raw
FITS cubes. These raw files are several hundred GB and are not included. The
compact PV caches needed by scripts 21 and 22 are included instead. To rebuild
the caches from raw data, pass the raw data paths explicitly:

```powershell
python 0_Data_preparation_2.py `
  --mwisp-fits <path-to-mosaic_12CO.fits> `
  --cfa-fits <path-to-CfA-FITS-cube> `
  --hi4pi-fits <path-to-HI4PI-CAR.fits>
```

The `MWISP_DATA_PATH` environment variable can also be used for the MWISP FITS
path.

### Local interactive 3D visualization

The local source file is built by:

```bash
cd code
python3 25_interactive_three_volume_viewer.py
```

Main inputs:

- `../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet`
- `../data/Bubbles.csv`
- `../results/superbubble_final_fit_parameters.csv`
- `../results/Open_superbubbles.csv` (optional simplified display table; the viewer falls back to the full table if this file is absent)

Outputs:

- `../results/interactive_3d_viewer.html`

The HTML embeds the uint8 dust-volume texture, open-superbubble geometry, and
bubble annotations. It loads Three.js 0.160.0 and its add-ons from `unpkg.com`,
so viewing requires network access. An optional `--output-bin <path>` exports a
separate copy of the volume payload; the HTML does not require that file.
In Code Ocean, the HTML is available after a run at
`/results/interactive_3d_viewer.html`.

### Output conventions

- `../data/` contains original catalogues and reusable non-generated caches.
- `../results/intermediate_output/` contains reproducible products, including
  generated 3D dust-map products, per-script CSV/JSON summaries, random trials,
  and diagnostic tables.
- `../results/` contains final manuscript-level products after a run.
- `../results/superbubble_final_fit_parameters.csv` is the
  authoritative 30-object geometry table generated by script 3.
- `../results/Open_superbubbles.csv` is the simplified catalogue table
  generated by script 20 after all upstream ingredients are available. Geometry
  and MCMC uncertainty columns come from script 3, shell and dust metrics from
  script 5, random-density p-values from script 8, dynamical ages and required
  SN counts from script 19, and the published evidence-class labels from the
  compact-catalogue class map encoded in script 20.
- `final_data_products_description.md` in this `code/` directory is a
  documentation file only; it explains the final data products and their
  columns and is not an analysis output.
- `pipeline_runtime_summary.png` in this `code/` directory is a static runtime
  summary figure for planning pipeline execution; it is
  not regenerated by the pipeline.
- `../results/SB_figures/` contains the six-panel maps for all 30
  objects after script 4 runs.
- `../results/figures/` contains manuscript, Extended Data, and
  diagnostic plots generated by scripts 3 and 5-24.

Most scripts set their working directory to `code/` and use code-relative paths.
The scripts use these direct relative paths both in Code Ocean and in the local capsule layout.


### Reproducibility caveats

Script 8 uses the full local raw 3D dust cube, including the cavity, shell,
and surroundings, to set the IAAFT Fourier amplitudes after missing-value
filling and optional block averaging. It preserves the density distribution
exactly and approximates the spectrum. The default is 8 iterations;
iteration count alone is not a convergence guarantee. See
[the method and diagnostic documentation in Script 8](8_OSBs_random_dust_test.py).

Each of the 30 targets also saves `8_SB{id}_original_vs_random.png` (PNG only)
in `../results/figures/random_density_diagnostics/`, using
the first realization and the standard 8-iteration setting. The comparison
has three rows of matched spatial slices (XY, XZ, YZ), with the original
on the left and the random realization on the right. Three panels beneath
them show the full-cube 3D power spectrum, power ratio and density distribution.
Figures retain concise titles, axes, units and legends; detailed method
descriptions and diagnostic numbers are recorded in the notes and metadata
for use in publication captions.
The six spatial slices use a shared `Spectral_r` colour scale and a display-only
Gaussian with sigma 1 voxel per in-plane axis on the actual test grid
(after optional rebinning), applied directly to each 2D slice.
Spectra, density distributions and saved arrays use
the cubes before display smoothing; details are in the method notes.
The arrays, spectral tables and metadata go to
`../results/intermediate_output/8_random_density_fluctuation_test/diagnostics/`.
Paths recorded in the diagnostic metadata are relative to `code/` for
portability; this path convention does not alter the scientific statistics.
From `code/`, run
`python 8_OSBs_random_dust_test.py --targets all --diagnostics-only --iaaft-max-iter 8`
to generate one comparison per target without rewriting Monte Carlo sample or
probability tables.
Use `--diagnostic-dir` and `--comparison-fig-dir` for separate diagnostic runs.

- `0_Data_preparation_2.py` processes the original HI4PI, CfA CO, and Galactic
  Plane Survey / MWISP FITS cubes. Those raw inputs are several hundred GB and
  are not included; the compact downstream PV caches are included.
- Script 8 is the longest standard analysis script because it runs 3D IAAFT
  surrogate-field tests.
- Script 9 automatically checks whether `dustmaps3d` is installed and, when it
  is available, locates the local 3D dust-map data package automatically. For
  normal execution, the standard `dustmaps3d` package cache is sufficient and
  no extra data-path argument is required.
- Scripts 14 and 15 can reuse cached orbit products by default. For a full orbit
  recomputation, run them with `--recompute-orbit --recompute-labels`; this can
  increase runtime.
- The composite final manuscript Fig. 1 is assembled from generated panels,
  as documented above. Script 23 writes manuscript Fig. 3 directly as one figure.
- If `../results/intermediate_output/` and `../results/` are deleted, the
  ordered script sequence can regenerate the reproducible products, provided the
  external `dustmaps3d` cache is available and the capsule `../data/`
  directory is kept.

### Companion online 3D dust-map platform

The companion online 3D dust-map platform for the open-superbubble
visualization is available at:

- https://nadc.china-vo.org/data/dustmaps/OSBs

### License

MIT.

### Contact

Prof. Haibo Yuan: yuanhb@bnu.edu.cn

Tao Wang: wt@mail.bnu.edu.cn

[返回顶部 / Back to top](#top) · [中文说明](#chinese)
