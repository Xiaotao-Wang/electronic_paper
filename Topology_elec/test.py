import numpy as np
import matplotlib.pyplot as plt
import mpmath

# 引入高精度计算：100位有效数字，确保 (C1/C2)^N 这种极端小/大数不丢精度
mpmath.mp.dps = 20

# 配置字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# =========================================================================
# 1. 严格对应你 LaTeX 里的参数
# =========================================================================
N_cells = 100  # 晶胞段数 N
C1 = mpmath.mpf('4e-11')  # 胞内弱耦合电容: 10 pF
C2 = mpmath.mpf('1e-9')  # 胞间强耦合电容: 1 nF
L = mpmath.mpf('8e-7')  # 节点接地电感

V_A0 = mpmath.mpf('1.0')  # 理想激励边界 V_A^0 = 1 V

# 比例系数 ratio = -C1 / C2
ratio = -C1 / C2

# =========================================================================
# 2. 直接套用你手算接好的有限长解析公式
# =========================================================================
# 初始化全网 200 个节点的复数电压数组
V_analytic = [mpmath.mpc('0', '0') for _ in range(2 * N_cells)]

for n in range(N_cells):
    # ---- 偶数点 A 节点 (索引 2n) ----
    # 公式: V_A^n = (-C1/C2)^n * V_A^0
    v_A = (ratio ** n) * V_A0
    V_analytic[2 * n] = v_A

    # ---- 奇数点 B 节点 (索引 2n+1) ----
    # 公式: V_B^n = (-C1/C2)^(2N - n + 1) * V_A^0
    # 注意：你公式里的末端项也可以统一写成指数 2N - n + 1 的形式
    v_B = (ratio ** (2 * N_cells - n + 1)) * V_A0
    V_analytic[2 * n + 1] = v_B

# =========================================================================
# 3. 抽取对数幅值供绘图
# =========================================================================
log_Vmag = []
for v in V_analytic:
    mag = mpmath.fabs(v)
    if mag < mpmath.mpf('1e-350'):
        log_Vmag.append(-350.0)
    else:
        log_Vmag.append(float(mpmath.log10(mag)))

# =========================================================================
# 4. 可视化渲染
# =========================================================================
plt.figure(figsize=(12, 6))
node_indices = np.arange(2 * N_cells)

# 绘制你接好的手算解析曲线
plt.plot(node_indices, log_Vmag, 'o-', color='tab:red', lw=1.2, ms=3, label='各个节点电压')

# 强行卡死纵坐标看清大跨度
plt.ylim(-400, 50)
tick_positions = np.arange(-350, 51, 50)
tick_labels = [f'$10^{{{i}}}$' for i in tick_positions]
plt.yticks(tick_positions, tick_labels)

plt.xlabel("一维链上点连续索引", fontsize=11)
plt.ylabel("对地电压幅值 |V| (V)", fontsize=11)
plt.title("C1=4e-11(F), C2=1e-9(F), L=8e-7(H), V_0=1(V), R_0=50(Ω)", fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11, loc='upper left')

plt.tight_layout()
plt.show()