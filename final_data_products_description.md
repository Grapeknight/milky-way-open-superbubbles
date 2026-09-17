# 最终数据产品说明 / Final Data Products Description

本文档说明可复现流程写入 `../results/` 的两份最终星表数据产品：

- `superbubble_final_fit_parameters.csv`：包含全部 30 个开放超泡的完整几何参数与来源记录，供分析脚本使用。
- `Open_superbubbles.csv`：经过精简与数值舍入的星表，便于论文展示和数据复用。运行可复现流程时，脚本 20 在上游几何、尘埃、随机密度、能量预算和超新星汇总产品生成后组装此表。

`superbubble_final_fit_parameters.csv` 是复现流程所采用的权威机器可读参数表；`Open_superbubbles.csv` 是面向科学结果展示的精简汇总表。

以下路径均以完整下载包的 `code/` 目录为基准；完整包保留 `code/data/results` 三目录结构。GitHub 仓库将 `code/` 的内容直接放在仓库根目录，完整数据与结果请从 [Releases](https://github.com/Grapeknight/milky-way-open-superbubbles/releases) 下载。

Paths below are relative to the `code/` directory of the complete download, which retains the `code/data/results` layout. The GitHub repository places the contents of `code/` directly at its root; download the complete data and results from [Releases](https://github.com/Grapeknight/milky-way-open-superbubbles/releases).

This documentation-only note describes the two catalogue-level final data
products written to `../results/` by the reproducible pipeline:

- `superbubble_final_fit_parameters.csv`: the full 30-object open-superbubble
  geometry and provenance table used by the analysis scripts.
- `Open_superbubbles.csv`: a simplified, rounded catalogue table intended for
  manuscript-facing use and convenient reuse. In a reproducible run, script 20
  assembles this compact table after the upstream geometry, dust, random-density,
  energy-budget, and SN-summary products are available.

`superbubble_final_fit_parameters.csv` is the authoritative machine-readable
table for reproducing the pipeline. `Open_superbubbles.csv` is the compact
science-facing summary.

## 坐标系与单位 / Coordinate System and Units

- 坐标采用以太阳为原点的银道笛卡尔坐标系，单位为 kpc。
- 坐标约定为 `x = d cos(b) cos(l)`、`y = d cos(b) sin(l)`、`z = d sin(b)`。
- 以 `_kpc`、`_pc`、`_deg` 结尾的字段，单位分别为 kpc、pc、度。
- MCMC 不确定度字段存储中位数相对于第 16 和第 84 百分位的单侧偏差。除非脚本明确另作说明，计算采用的仍是非 MCMC 的最终几何参数字段。

- Coordinates are heliocentric Galactic Cartesian coordinates in kpc.
- The convention is:
  - `x = d cos(b) cos(l)`
  - `y = d cos(b) sin(l)`
  - `z = d sin(b)`
- Columns ending in `_kpc` are in kpc; columns ending in `_pc` are in pc; columns
  ending in `_deg` are in degrees.
- MCMC uncertainty columns store one-sided 16th/84th percentile offsets around
  the median. The adopted values used by the scripts remain the non-MCMC final
  geometry columns unless a script explicitly states otherwise.

## 目标类型与几何分类 / Object and Geometry Classes

### 完整表中的 `mark` / `mark` in `superbubble_final_fit_parameters.csv`

| 取值 / Value | 含义 / Meaning | 几何约定 / Geometry convention |
|---:|---|---|
| `1` | 北向开放超泡<br>North-open OSB | 按照脚本所用的椭球帽约定，开口帽沿局部及全局 Z 的负方向。<br>Open cap points toward negative local/global Z in the cap convention used by the scripts. |
| `2` | 南向开放超泡<br>South-open OSB | 按照脚本所用的椭球帽约定，开口帽沿局部及全局 Z 的正方向。<br>Open cap points toward positive local/global Z in the cap convention used by the scripts. |
| `3` | 穿盘超泡<br>Disk-penetrating OSB | 通常采用椭圆柱模型。<br>Usually represented by an elliptical-cylinder model. |

### 完整表中的 `shape` / `shape` in `superbubble_final_fit_parameters.csv`

| 取值 / Value | 含义 / Meaning |
|---|---|
| `ellipsoid` | 椭球目标。表中给出拟合椭球，以及 `xy_plane_*` 字段中的尘埃平面帽底参数。<br>Ellipsoidal target. The table gives the fitted ellipsoid and the dust-plane cap-base parameters in `xy_plane_*`. |
| `cylinder` | 椭圆柱目标。完整表将 `c_radius_kpc` 设为 `0`，以标识柱体约定；各脚本自行定义所需的垂直半高。<br>Elliptical-cylinder target. `c_radius_kpc` is set to `0` in the full table to mark the cylinder convention; individual scripts define the vertical half-height they need. |

### 精简表中的 `type` / `type` in `Open_superbubbles.csv`

| 取值 / Value | 含义 / Meaning |
|---|---|
| `N` | 北向开放超泡，对应 `mark = 1`。<br>North-open OSB, corresponding to `mark = 1`. |
| `S` | 南向开放超泡，对应 `mark = 2`。<br>South-open OSB, corresponding to `mark = 2`. |
| `P` | 穿盘超泡，对应 `mark = 3`。<br>Disk-penetrating OSB, corresponding to `mark = 3`. |

### 精简表中的 `class` / `class` in `Open_superbubbles.csv`

`class` 是精简星表使用的证据分类。完整表保留了组装精简星表所依据的数值诊断量。

`class` is the compact evidence class used for the simplified catalogue. The
full table retains the numerical diagnostics from which the simplified
catalogue was assembled.

## 完整表 / Full Table: `superbubble_final_fit_parameters.csv`

30 行，69 列。 / Rows: 30. Columns: 69.

| 字段 / Column | 说明 / Description |
|---|---|
| `id` | 整套数据与代码中使用的开放超泡（OSB）编号。<br>OSB identifier used throughout the package. |
| `version` | 最终参数表的版本标签。<br>Final-parameter table version label. |
| `mark` | 开口或形态类别，定义见上文。<br>Opening or morphology class, defined above. |
| `shape` | 三维几何模型，定义见上文。<br>Three-dimensional geometry model, defined above. |
| `xy_model` | XY 截面采用椭圆还是圆化模型。<br>Whether the XY cross-section is elliptical or circularized. |
| `circular_xy` | 最终拟合是否采用圆化 XY 截面的布尔标记。<br>Boolean flag for the circularized-XY final fit. |
| `potential_circular_xy_by_axis_ratio` | 初始轴比规则是否建议圆化。<br>Whether the initial axis-ratio rule suggested circularization. |
| `final_parameter_estimator` | 所采用最终参数的来源，例如 `least_squares` 或 `manual_override`。<br>Source of the adopted final parameters, such as `least_squares` or `manual_override`. |
| `final_estimator_reason` | 采用该最终参数来源的简短原因记录。<br>Short provenance note for why the final parameter source was adopted. |
| `initial_free_ab_axis_ratio` | 初始、未约束的自由 XY 解的轴比。<br>Axis ratio from the unconstrained initial/free XY solution. |
| `fit_ab_axis_ratio` | 采用的最终 XY 轴比 `a/b`；圆化目标取值为 `1`。<br>Adopted final XY axis ratio `a/b`; circularized targets have value `1`. |
| `least_squares_a_kpc` | 展示精简前的最小二乘半轴 `a`。<br>Least-squares semi-axis `a` before any reporting simplification. |
| `least_squares_b_kpc` | 展示精简前的最小二乘半轴 `b`。<br>Least-squares semi-axis `b` before any reporting simplification. |
| `least_squares_c_kpc` | 最小二乘半轴 `c`；柱体目标采用上文定义的柱体约定。<br>Least-squares semi-axis `c`; cylinder targets use the cylinder convention described above. |
| `least_squares_angle_deg` | 最小二乘 XY 位置角。<br>Least-squares XY position angle. |
| `center_x_kpc` | 采用的最终中心 X 坐标。<br>Adopted final center X coordinate. |
| `center_y_kpc` | 采用的最终中心 Y 坐标。<br>Adopted final center Y coordinate. |
| `center_z_kpc` | 采用的最终中心 Z 坐标。<br>Adopted final center Z coordinate. |
| `a_radius_kpc` | 采用的最终局部半轴 `a`。<br>Adopted final local `a` semi-axis. |
| `b_radius_kpc` | 采用的最终局部半轴 `b`。<br>Adopted final local `b` semi-axis. |
| `c_radius_kpc` | 采用的最终局部椭球半轴 `c`；取值 `0` 表示采用柱体约定。<br>Adopted final local `c` semi-axis for ellipsoids; `0` marks the cylinder convention. |
| `angle_deg` | 局部 `a` 轴的位置角，从全局 `+X` 方向逆时针测量。<br>Position angle of the local `a` axis, measured counter-clockwise from global `+X`. |
| `fit_std_pc` | 所选首次接触分子云相对拟合壳层的残差标准差。<br>Standard deviation of selected first-contact molecular-cloud residuals from the fitted shell. |
| `C_local` | 局部一致性统计量；人工处理的目标可能为空。<br>Local consistency statistic; may be blank for manually handled targets. |
| `catalog_cloud_count` | 源星表中的分子云数量。<br>Number of molecular clouds in the source catalogue. |
| `local_bubble_removed_cloud_count` | 剔除本地泡边界分子云后的数量。<br>Number after removing Local Bubble boundary clouds. |
| `candidate_cloud_count` | 预选候选分子云数量。<br>Number of preselected candidate clouds. |
| `selected_cloud_count` | 用于拟合的首次接触分子云数量。<br>Number of first-contact clouds used in the fit. |
| `ray_count` | 首次接触搜索使用的射线数量。<br>Number of rays used by the first-contact search. |
| `ray_hit_count` | 保留的射线命中记录数量。<br>Number of retained ray-hit records. |
| `hit_ray_count` | 至少命中一次的射线数量。<br>Number of rays with at least one hit. |
| `catalog_filter` | 分子云星表的筛选规则。<br>Molecular-cloud catalogue filtering rule. |
| `xy_search_scale` | 无量纲 XY 搜索尺度控制参数。<br>Dimensionless XY search-scale control. |
| `fit_center_bounds_scale` | 拟合时使用的无量纲中心边界尺度。<br>Dimensionless center-bound scale used during fitting. |
| `search_scale` | 首次接触阶段使用的无量纲搜索尺度控制参数。<br>Dimensionless search-scale control used by the first-contact stage. |
| `radius_tolerance_pc` | 首次接触径向容差。<br>First-contact radial tolerance. |
| `ray_cone_aperture_deg` | 用于将分子云关联到射线的射线锥孔径。<br>Ray-cone aperture used to associate clouds with a ray. |
| `max_hits_per_ray` | 每条射线最多保留的命中次数。<br>Maximum number of retained hits per ray. |
| `fit_cloud_seq_list` | 拟合使用的分子云序列编号。<br>Molecular-cloud sequence IDs used by the fit. |
| `fit_clouds_json` | JSON 编码的所选分子云记录及局部坐标系拟合诊断信息。<br>JSON-encoded selected cloud records and local-frame fit diagnostics. |
| `xy_plane_center_x_kpc` | 代表性椭球帽底面或柱体底面椭圆的 X 坐标。<br>X coordinate of the representative cap-base or cylinder-base ellipse. |
| `xy_plane_center_y_kpc` | 代表性椭球帽底面或柱体底面椭圆的 Y 坐标。<br>Y coordinate of the representative cap-base or cylinder-base ellipse. |
| `xy_plane_z_kpc` | 代表性尘埃平面椭圆的 Z 高度。<br>Z height of the representative dust-plane ellipse. |
| `xy_plane_a_kpc` | 代表性尘埃平面椭圆的半轴 `a`。<br>Semi-axis `a` of the representative dust-plane ellipse. |
| `xy_plane_b_kpc` | 代表性尘埃平面椭圆的半轴 `b`。<br>Semi-axis `b` of the representative dust-plane ellipse. |
| `xy_plane_angle_deg` | 代表性尘埃平面椭圆的位置角。<br>Position angle of the representative dust-plane ellipse. |
| `xy_plane_definition` | 代表性椭圆的来源，例如 `cap_base_ellipse` 或 `cylinder_base_ellipse`。<br>Provenance of the representative ellipse, such as `cap_base_ellipse` or `cylinder_base_ellipse`. |
| `mcmc_center_x_median_kpc` | 中心 X 坐标的 MCMC 中位数。<br>MCMC median for center X. |
| `mcmc_center_x_err_minus_kpc` | 中心 X 坐标的下侧单侧不确定度。<br>Lower one-sided uncertainty for center X. |
| `mcmc_center_x_err_plus_kpc` | 中心 X 坐标的上侧单侧不确定度。<br>Upper one-sided uncertainty for center X. |
| `mcmc_center_y_median_kpc` | 中心 Y 坐标的 MCMC 中位数。<br>MCMC median for center Y. |
| `mcmc_center_y_err_minus_kpc` | 中心 Y 坐标的下侧单侧不确定度。<br>Lower one-sided uncertainty for center Y. |
| `mcmc_center_y_err_plus_kpc` | 中心 Y 坐标的上侧单侧不确定度。<br>Upper one-sided uncertainty for center Y. |
| `mcmc_center_z_median_kpc` | 中心 Z 坐标的 MCMC 中位数。<br>MCMC median for center Z. |
| `mcmc_center_z_err_minus_kpc` | 中心 Z 坐标的下侧单侧不确定度。<br>Lower one-sided uncertainty for center Z. |
| `mcmc_center_z_err_plus_kpc` | 中心 Z 坐标的上侧单侧不确定度。<br>Upper one-sided uncertainty for center Z. |
| `mcmc_a_median_kpc` | 半轴 `a` 的 MCMC 中位数。<br>MCMC median for semi-axis `a`. |
| `mcmc_a_err_minus_kpc` | 半轴 `a` 的下侧单侧不确定度。<br>Lower one-sided uncertainty for semi-axis `a`. |
| `mcmc_a_err_plus_kpc` | 半轴 `a` 的上侧单侧不确定度。<br>Upper one-sided uncertainty for semi-axis `a`. |
| `mcmc_b_median_kpc` | 半轴 `b` 的 MCMC 中位数。<br>MCMC median for semi-axis `b`. |
| `mcmc_b_err_minus_kpc` | 半轴 `b` 的下侧单侧不确定度。<br>Lower one-sided uncertainty for semi-axis `b`. |
| `mcmc_b_err_plus_kpc` | 半轴 `b` 的上侧单侧不确定度。<br>Upper one-sided uncertainty for semi-axis `b`. |
| `mcmc_c_median_kpc` | 半轴 `c` 或展示所用垂直尺度的 MCMC 中位数。<br>MCMC median for semi-axis `c` or the reporting vertical scale. |
| `mcmc_c_err_minus_kpc` | `c` 的下侧单侧不确定度。<br>Lower one-sided uncertainty for `c`. |
| `mcmc_c_err_plus_kpc` | `c` 的上侧单侧不确定度。<br>Upper one-sided uncertainty for `c`. |
| `mcmc_angle_median_deg` | 位置角的 MCMC 中位数。<br>MCMC median for the position angle. |
| `mcmc_angle_err_minus_deg` | 位置角的下侧单侧不确定度。<br>Lower one-sided uncertainty for the position angle. |
| `mcmc_angle_err_plus_deg` | 位置角的上侧单侧不确定度。<br>Upper one-sided uncertainty for the position angle. |
| `mcmc_acceptance_fraction_mean` | 该目标的平均 MCMC 接受率。<br>Mean MCMC acceptance fraction for the target. |

## 精简表 / Simplified Table: `Open_superbubbles.csv`

30 行，41 列。 / Rows: 30. Columns: 41.

此表用于简明报告科学结果，数值经过舍入以便阅读。脚本 20 从当前分析产品重新生成该表：几何参数与 MCMC 不确定度来自脚本 3；壳层和尘埃指标来自脚本 5；尘埃联合 p 值来自脚本 8；动力学年龄与所需超新星数目来自脚本 19；精简证据分类标签来自脚本 20 内的分类映射。无需预先提供静态的 `Open_superbubbles.csv` 输入文件。本表不包含内部射线与分子云来源记录，也不包含 MCMC 运行记录。

This table is a compact reporting catalogue. Values are rounded for readability.
Script 20 regenerates it from the current analysis products: geometry and MCMC
uncertainties from script 3, shell and dust metrics from script 5, joint dust
p-values from script 8, dynamical ages and required supernova counts from script
19, and the compact evidence-class labels from the class map encoded in script
20. No static `Open_superbubbles.csv` input is required.
The table omits internal ray/cloud provenance and MCMC bookkeeping.

| 字段 / Column | 说明 / Description |
|---|---|
| `id` | 整套数据与代码中使用的开放超泡（OSB）编号。<br>OSB identifier. |
| `x_c_kpc` | 报告的中心 X 坐标。<br>Reported center X coordinate. |
| `x_c_err_minus_kpc` | 中心 X 坐标的下侧不确定度。<br>Lower uncertainty for center X. |
| `x_c_err_plus_kpc` | 中心 X 坐标的上侧不确定度。<br>Upper uncertainty for center X. |
| `y_c_kpc` | 报告的中心 Y 坐标。<br>Reported center Y coordinate. |
| `y_c_err_minus_kpc` | 中心 Y 坐标的下侧不确定度。<br>Lower uncertainty for center Y. |
| `y_c_err_plus_kpc` | 中心 Y 坐标的上侧不确定度。<br>Upper uncertainty for center Y. |
| `z_c_kpc` | 报告的中心 Z 坐标。对于穿盘目标，此值遵循精简报告约定，不采用 `c_radius_kpc` 中的柱体标记。<br>Reported center Z coordinate. For disk-penetrating targets, this follows the simplified reporting convention rather than the cylinder marker in `c_radius_kpc`. |
| `z_c_err_minus_kpc` | 中心 Z 坐标的下侧不确定度。<br>Lower uncertainty for center Z. |
| `z_c_err_plus_kpc` | 中心 Z 坐标的上侧不确定度。<br>Upper uncertainty for center Z. |
| `a_kpc` | 报告的半轴 `a`。<br>Reported semi-axis `a`. |
| `a_err_minus_kpc` | `a` 的下侧不确定度。<br>Lower uncertainty for `a`. |
| `a_err_plus_kpc` | `a` 的上侧不确定度。<br>Upper uncertainty for `a`. |
| `b_kpc` | 报告的半轴 `b`。<br>Reported semi-axis `b`. |
| `b_err_minus_kpc` | `b` 的下侧不确定度。<br>Lower uncertainty for `b`. |
| `b_err_plus_kpc` | `b` 的上侧不确定度。<br>Upper uncertainty for `b`. |
| `c_kpc` | 报告的垂直半轴或展示所用的垂直尺度。<br>Reported vertical semi-axis or reporting vertical scale. |
| `c_err_minus_kpc` | `c` 的下侧不确定度。<br>Lower uncertainty for `c`. |
| `c_err_plus_kpc` | `c` 的上侧不确定度。<br>Upper uncertainty for `c`. |
| `pa_deg` | 报告的位置角。<br>Reported position angle. |
| `pa_err_minus_deg` | 位置角的下侧不确定度。<br>Lower uncertainty for position angle. |
| `pa_err_plus_deg` | 位置角的上侧不确定度。<br>Upper uncertainty for position angle. |
| `class` | 精简证据分类。<br>Compact evidence class. |
| `type` | 精简形态类别：`N`、`S` 或 `P`，定义见上文。<br>Compact morphology type: `N`, `S`, or `P`, as defined above. |
| `a_compromise_kpc` | 报告及派生计算使用的折衷投影半轴 `a`。<br>Projected compromise semi-axis `a` used by reporting/derived calculations. |
| `b_compromise_kpc` | 报告及派生计算使用的折衷投影半轴 `b`。<br>Projected compromise semi-axis `b` used by reporting/derived calculations. |
| `r_eff_compromise_kpc` | 折衷投影的有效半径。<br>Effective projected compromise radius. |
| `pa_compromise_deg` | 折衷投影椭圆采用的位置角。<br>Position angle used with the compromise projected ellipse. |
| `fit_std_pc` | 所选首次接触分子云相对拟合壳层的残差标准差。<br>First-contact shell fit residual scatter. |
| `shell_radius_mean_abs_offset_over_R_eq` | 以等效半径归一化的壳层半径平均绝对偏移。<br>Mean absolute shell-radius offset normalized by the equivalent radius. |
| `shell_ridge_uPeak_curvature_rms` | 归一化半径下壳层脊峰位置的曲率均方根统计量。<br>RMS curvature statistic of the shell-ridge peak in normalized radius. |
| `dust_shell_ridge_peak_inner_ratio` | 尘埃壳层脊峰与内部空腔尘埃水平之比。<br>Ratio of dust-shell ridge peak to inner-cavity dust level. |
| `dust_inner_mean_mag_kpc` | 内部空腔尘埃密度均值，单位为 mag/kpc。<br>Mean inner-cavity dust density in mag/kpc. |
| `dust_inner_std_mag_kpc` | 内部空腔尘埃密度标准差，单位为 mag/kpc。<br>Standard deviation of inner-cavity dust density in mag/kpc. |
| `p_dust_joint3` | 精简证据汇总采用的尘埃联合显著性 p 值。<br>Joint dust-significance p-value used for the compact evidence summary. |
| `age_dyn_myr` | 动力学年龄估计，单位为 Myr。<br>Dynamical age estimate in Myr. |
| `age_myr_err_minus` | 动力学年龄的下侧不确定度。<br>Lower uncertainty for the dynamical age. |
| `age_myr_err_plus` | 动力学年龄的上侧不确定度。<br>Upper uncertainty for the dynamical age. |
| `n_sn_required` | 估计的所需超新星数目。<br>Estimated number of required supernovae. |
| `n_sn_required_err_minus` | 所需超新星数目的下侧不确定度。<br>Lower uncertainty for the required supernova count. |
| `n_sn_required_err_plus` | 所需超新星数目的上侧不确定度。<br>Upper uncertainty for the required supernova count. |

## 几何重建说明 / Geometry Reconstruction Notes

对于椭球，首先将点旋转到局部坐标系：

For an ellipsoid, first rotate a point into the local frame:

```python
import numpy as np

def rotate_to_local(points, center, angle_deg):
    q = np.asarray(points, dtype=float) - np.asarray(center, dtype=float)
    th = np.deg2rad(-angle_deg)
    c, s = np.cos(th), np.sin(th)
    rot = np.array([[c, -s, 0.0],
                    [s,  c, 0.0],
                    [0.0, 0.0, 1.0]])
    return q @ rot.T
```

局部 `x'` 轴对应半轴 `a`，局部 `y'` 轴对应半轴 `b`，局部 `z'` 轴与全局 Z 轴平行。完整椭球表面满足下式：

The local `x'` axis corresponds to `a`, the local `y'` axis to `b`, and the
local `z'` axis is parallel to global Z. The full ellipsoid surface is:

```text
u_full = sqrt((x'/a)^2 + (y'/b)^2 + (z'/c)^2) = 1
```

在尘埃壳层统计中，椭球类开放超泡使用由 `xy_plane_*` 定义的自相似椭球帽几何。典型壳层与内部掩膜的定义见执行相应统计的各个脚本。

For dust-shell statistics, ellipsoidal OSBs use the self-similar cap geometry
defined by `xy_plane_*`. Typical shell/interior masks are documented in the
individual scripts that compute the corresponding statistics.
