# 任意形状的外部再电离历史

本分支在上游 CAMB **2.0.5**（基准提交
`1d1d3d6d09e499291651855d4d8286419030a32b`）上增加
`camb.reionization.TabulatedReionization`，包版本为 **2.0.5+tabulated.1**。
它接受表格或一般 Python 函数描述的、空间均匀的复合后电子历史，不要求历史全局单调，
也没有人为设置的 `z=30` 或 `z=50` 上限。

**接口能够表示某条历史，不等于这条历史有自洽的物理来源，或已证明其数值收敛。**
本说明对应当前实现；背景和 CMB 谱的实际验收数值、所用精度设置与剩余限制，
由本分支验证报告另行记录，不以本文示例代替验证结果。

## 适用范围与科学边界

- 用途：给定外部电子历史，计算其均匀背景 Thomson 散射和线性 CMB 响应。
  可表示高红移部分电离、早期电离后再复合、多段或非单调的有效历史。
- 不求解电离源、复合率、辐射输运、能量沉积或气体加热；不自动产生自洽的温度和 21-cm 历史。
  也不表示空间斑块性、patchy kSZ 或再电离反馈。
- 不替换复合模块。越靠近复合时期，越不能仅凭“红移未被拒绝”判断模型适用；
  若修改了复合或温度演化，需另外审查该模块和外部物理模型的耦合。
- CAMB v2 默认让完整电子历史进入 photon opacity，但不把该历史自动当成自洽加热史。
  `include_heating` 默认关闭；启用后也只是粗略的均匀加热诊断，不能据其推断真实 IGM 热史或小尺度物质反馈。

## 先选对物理定义

CAMB 使用的电子丰度是电子数与氢核数之比，不是氢电离分数：

\[
x_e(z)\equiv n_e(z)/n_{\rm H}(z),\qquad
f_{\rm He}\equiv n_{\rm He}/n_{\rm H}
=\frac{Y_{\rm He}}{(m_{\rm He}/m_{\rm H})(1-Y_{\rm He})}.
\]

实现使用 CAMB 的原子质量比常数，而不是将质量比精确写成 4。
记 `r(z)` 为当前宇宙学参数和复合模块给出的残余电子丰度；
`F=1+f_He` 为氢完全电离、氦单电离的丰度。
以下公式先不含单独的 He II → He III 项，也暂不考虑表格外的连接区。

| `history_mode` | 输入值含义 | 实际送入 CAMB 的电子丰度 |
| --- | --- | --- |
| `ionization`（默认） | `t(z)`：主再电离组分，以每氢核电子数计，不含未电离区域的 residual | `q=t/F`，`x_e=(1-q)r+qF=t+(1-t/F)r` |
| `ionized_fraction` | `q(z)`：主再电离组分的归一化有效比例 | `x_e=(1-q)r+qF` |
| `total_xe` | `t(z)`：已经含 residual 的总电子丰度 | `x_e=t`；不再添加 residual |

因此 `ionization` **不是**简单的 `x_e=r+t`，`ionized_fraction` **不是**总 `x_e`。
前两种模式隐含氢与氦第一次电离共享同一个有效比例；若需要独立的氦历史，
请在外部构造正确的总丰度，使用 `total_xe`。

输入需为等长、有限、实数、非负的一维数组，至少两个不同红移；输入会按红移排序，
不允许重复红移或负红移。`ionized_fraction` 限于 `[0,1]`。
其余模式的电子预算在 CAMB 初始化时依据实际 `YHe` 检查：
`ionization` 上限为 `F`；`total_xe` 在不另加 HeIII 时上限为 `1+2f_He`，
另加 HeIII 时输入上限为 `F`。不把不合法值静默裁剪到范围内。

### 氦只加一次

`include_helium_fullreion` 控制 CAMB 单独的低红移 He II → He III 步骤。
前两种模式默认开启，`total_xe` 默认关闭。
开启时，最终结果还会加上该步骤的电子贡献 `h(z)`；表外连接到的是 `r(z)+h(z)`。
若输入已包含 HeIII，必须关闭它，避免重复计数。
这里的 `helium_redshift` 等参数仍控制原有 CAMB 氦步骤；本接口并未独立求解氦的原子过程。
自定义氦参数时须有限、宽度为正，且开始红移须大于
`max(0, helium_redshift-8*helium_delta_redshift)`；无效源区间会在初始化时报错。

### 从非均匀模拟转换时注意平均方式

均匀 opacity 严格需要与平均电子密度一致的
`<n_e>/<n_H>`；一般不等于局域电离分数的简单体积平均。
直接把模拟的体积平均 `1-<x_HI>` 当成 `q`，隐含忽略电离状态与密度关联的近似，
还需核对模拟的 residual 和氦约定。该转换不由此接口自动保证。

## 两端和复合历史如何衔接

设原始表覆盖 `[z_min,z_max]`，高红移连接宽度为 `D=transition_width`。

- **低红移端**：输入值保持为第一个节点值，不强制完全电离。
  这表示 `q` 或 `t` 的常数延伸；前两种模式的实际 `x_e` 仍通过公式含有随红移变化的 `r(z)`，
  可选氦项也仍正常作用。
- **表内**：保留原始节点值，不因为需要连接复合史而改动表内节点。
- **高红移连接区** `z_max < z < z_max+D`：令
  `u=(z-z_max)/D`、`w=10u^3-15u^4+6u^5`。
  前两种模式将末端有效比例替换为 `q_end*(1-w)`，再代入上面的 residual 混合公式；
  `total_xe` 使用 `x_e=(1-w)t_end+w*r(z)`。
  最后按配置添加一次 `h(z)`。
- **连接区以上**：使用 `r(z)` 加可选 `h(z)`，不将末端非零电离永久延伸到任意高红移。

宽度规则：

1. `ionization` 或 `ionized_fraction` 的末端输入恰为零时，省略宽度等价于 `D=0`，直接接上 residual。
2. 末端非零但不大于 `boundary_tolerance`（默认 `1e-6`）时，省略宽度会自动选
   `D=z[-1]-z[-2]`。这是追加在表外的一段连接，不会把原节点截成零。
   容差检查的是**该模式的原始输入单位**，不是统一换算后的总电子丰度。
3. 末端超过该容差时，必须显式给正宽度，或先延长输入历史到足够小的尾部；否则报错。
4. `total_xe` **始终要求显式正宽度**，即使末端数值很小。因为它含 residual，
   不能把“小丰度”自动解释成“已经接上当前宇宙学的复合解”。

选择正宽度相当于明确增加了一段模型历史，会改变光深和谱；它不是免费的数值修饰。
对于物理来源已知的尾部，优先在表内完整提供，而不是用连接区替代未知物理。

## 插值：受限的 C2 Bernstein，不是 PCHIP

每两个节点间使用五次 Bernstein 多项式。Python 先从 clamped cubic spline
估计共享的一阶、二阶节点导数，再限制导数，使各段 Bernstein 控制点按该段端点方向排列。
Fortran 使用 de Casteljau 算法求值。

这个实现的目的有两个：共享导数使曲线为 **C2**；每段都限制在两端数据范围内，避免过冲。
两端的一阶、二阶导数设为零，与常数输入延伸衔接。历史可以整体非单调；
保持的是每段的局部方向和端点范围。因此 `[0,1]` 的 `q` 节点之间也不会凭插值跑出该范围。

**这不是 SciPy PCHIP，也不是文献 FlexKnot 的原样实现。** PCHIP 通常为 C1；
这里为了热力学插值相关导数的平滑性选用受限 C2 表示。
五次多项式不保证全局五阶精度；导数限制、稀疏节点和尖锐特征都会影响误差。
表格足够稠密与否，必须通过加密节点后比较实际背景、光深和谱来判断。

## 示例一：z≈55 的非单调电离峰

这是演示接口的假设历史，不是观测最佳拟合模型或完整的早期电离物理模型。
表延伸至 `z=80` 且末端为零；因此不需人为指定高红移连接宽度。

```python
import numpy as np
import camb
from camb.reionization import TabulatedReionization

z = np.array([0., 5., 7., 9., 15., 30., 45., 50., 55., 60., 65., 75., 80.])
q = np.array([1., 1., .65, .12, .01, .002, .002, .01, .03, .01, .002, 0., 0.])
options = dict(history_mode="ionized_fraction", include_helium_fullreion=True)

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.12, mnu=0.06, YHe=0.245)
pars.InitPower.set_params(As=2.1e-9, ns=0.965)
pars.Reion = TabulatedReionization().set_table(z, q, **options)

# A short background check, not a spectrum convergence test.
background = camb.get_background(pars)
z_check = np.array([0., 8., 30., 55., 80., 85.])
history = background.get_background_redshift_evolution(z_check, ["x_e"])
print(camb.__version__, camb.__file__)
print(np.column_stack((z_check, history["x_e"])))

# Store the actual input, not a pickled Fortran object.
np.savez("example_history.npz", redshifts=z, values=q)
```

请另外保存 `options`、完整宇宙学/复合参数和版本提交号；重载时从数组重新 `set_table`。
一般历史没有唯一的 `zrei`，`get_zre`、`set_zrei` 和 `set_tau` 会明确报错。
`tau` 应由历史导出，而不是再给定一个彼此不一致的自由输入。

## 示例二：解析函数输入

`set_function` 先在 Python 自适应采样，再使用与表格接口完全相同的 Fortran 表示；
没有在积分热循环里反复执行 Python callback。

```python
import numpy as np
from camb.reionization import TabulatedReionization

def q_of_z(z):
    return 0.5 * (1.0 + np.tanh((8.0 - z) / 0.8))

reion = TabulatedReionization()
diagnostics = reion.set_function(
    q_of_z, 0.0, 15.0,
    history_mode="ionized_fraction",
    transition_width=2.0,
    atol=1e-6, rtol=1e-4, max_dz=0.25,
    breakpoints=[7.0, 8.0, 9.0], max_nodes=10000,
)
print(diagnostics)
```

这里显式声明 `transition_width`；若省略，则函数在 `z_max` 的返回值必须满足上述小尾部规则。
诊断的含义与默认值为：

- `atol=1e-6`、`rtol=1e-4`：各区间的 1/4、1/2、3/4 位置比较真实函数与插值，
  检查 `abs(error) <= atol + rtol*abs(function(z))`；不是导数误差或 CMB 谱误差。
- `max_dz=1.0`：初始网格最大间隔。需按已知物理特征尺度缩小，而不是只依赖误差探针。
- `breakpoints=()`：将已知峰、过渡等位置加入节点；它不会自动推断特征宽度，仍须用间隔解析其宽度。
- `max_nodes=10000`：节点预算；预算耗尽而未满足探针条件会报错，不静默接受。
- 返回 `n_nodes`、`max_probe_error`、容差、最大初始间隔、迭代次数，以及历史模式和实际连接宽度。

任何有限采样都可能漏过黑箱函数中未告知的极窄峰；“探针全部通过”不是任意函数的一致误差证明。
函数须确定性、有限、非负，可接受数组输入或逐点标量输入。非平坦端点与这里的零导数端点约定不同，
可能需要额外加密。

需要保存函数得到的**实际采样数组**时，使用公共 `get_table()` 接口：

```python
table = reion.get_table()
np.savez("sampled_function_history.npz", **table)
# table contains independent copies of redshifts and x_e, in the INPUT mode's units.
```

保存数组、全部 `set_table` 选项及模型参数，
比保存一个匿名函数或只存采样容差更可复现。
本类不支持 pickle，也不能用 CAMB `.ini` 序列化其历史数组；这两条路径会报错。

## 安装隔离与数值验收

不要覆盖现有生产 `cosmo` 环境或原项目的 patched CAMB。应在新环境、新进程中构建本工作树，
并检查 Python 和 Fortran 版本以及实际导入路径。例如，具备 Python 3.10+ 和 gfortran 后：

```sh
git clone --recurse-submodules --branch codex/tabulated-reionization-v2 https://github.com/SeanZhang226/CAMB.git CAMB-tabulated
python3 -m venv ~/.venvs/camb-tabulated-v2
~/.venvs/camb-tabulated-v2/bin/python -m pip install -e ./CAMB-tabulated
~/.venvs/camb-tabulated-v2/bin/python -c 'import camb; print(camb.__version__, camb.__file__)'
```

首次安装会构建 Fortran。修改 Fortran 后，在仓库目录用该环境的
`python setup.py make` 重建；不要在已载入旧共享库的 Python 进程中继续验证新代码。
上述路径为本地示例，其他机器替换为自己的 checkout。

验收至少分开检查：

1. 输入语义、氦不重复、电子预算和两端规则。
2. 实际 CAMB `x_e(z)` 与目标曲线，包括峰、连接区和 residual，而不只看 Python 插值。
3. 用一致单位独立积分光深；采用物理氢密度时，
   `tau = c*sigma_T*n_H0*integral[(1+z)^2*x_e(z)/H(z) dz]`。
   应报告积分区间，并区分总电子光深与相对 residual 的额外贡献。
4. 分别加密输入表、背景网格和源函数时间采样。`min_steps`（默认 200）及
   `timestep_boost` 只约束再电离源函数区间的最少采样，不保证热力学网格或最终谱精度；
   还需检查 `BackgroundTimeStepBoost`、`AccuracyBoost` 等控制的收敛。
   本分支还在每个有变化的表格区间增加 `source_samples_per_interval` 个源采样（默认 16），
   并解析显式连接区及可选氦步骤；实际采样数也随 `AccuracyBoost*TimeStepBoost` 增加。
   这避免只用全局再电离起止区间而漏过局部峰，但不等于任意窄峰都被背景网格解析。
   `boundary_tolerance` 只决定末端连接规则，不用来丢弃小振幅的表内结构。
5. 比较相同光深但高红移分布不同的历史，检查 EE 至至少 `ell=100` 及相关 TT/TE；
   不能仅凭 tau 一致判定科学响应一致。TE 过零处采用绝对误差或按 `sqrt(TT*EE)` 归一化。

本文不宣称示例的节点数、采样默认值或任何单一红移范围已经达到给定谱精度。
实际背景/谱验收、精度设置及原生指数模型对照的限制见[验证报告](tabulated_validation.md)。

## 参考资料

- Antony Lewis, *CAMB v2: cosmological power spectra for high-precision surveys*,
  JCAP 10 (2026) 013：[arXiv:2607.14854](https://arxiv.org/abs/2607.14854)。
  参见 V.2、V.3 和附录 B/C，特别是复合、电离定义及非自洽均匀加热的限制。
- Chen Heinrich and Wayne Hu, *Does Planck 2015 polarization favor high redshift reionization?*,
  Phys. Rev. D 98, 063514 (2018)：[arXiv:1802.00791](https://arxiv.org/abs/1802.00791)。
  包含 `z_max=30,40,50` 的历史和极化形状响应；其数据时期和结论不能等同当前观测证据。
- Marius Millea and François Bouchet, *Cosmic Microwave Background Constraints in Light of Priors Over Reionization Histories*,
  A&A 617, A96 (2018)：[arXiv:1804.08476](https://arxiv.org/abs/1804.08476)。
  讨论历史先验、早期电离推断及使用 PCHIP 的 FlexKnot；本接口不等同该统计模型。
- Hongwan Liu, Gregory W. Ridgway and Tracy R. Slatyer,
  *DarkHistory: A code package for calculating modified cosmic ionization and thermal histories with dark matter and other exotic energy injections*,
  Phys. Rev. D 101, 023530 (2020)：[arXiv:1904.09296](https://arxiv.org/abs/1904.09296)。
  展示具体能量注入模型中同时演化电离和温度的重要性。
- [SciPy PchipInterpolator 官方说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html)：
  仅用于对照 C1、保持形状的 PCHIP；不是本分支 Bernstein 算法的误差定理。
