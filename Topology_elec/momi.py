import numpy as np
import matplotlib.pyplot as plt
import mpmath

# 引入高精度计算：手动指定保留 100 位有效数字，突破双精度浮点数机器极限
mpmath.mp.dps = 100

# 配置字体
plt.rcParams['font.sans-serif'] = ['SimHei']     # 配置中文字体（黑体）
plt.rcParams['axes.unicode_minus'] = False        # 3. 彻底把 Unicode 负号开关关掉，用普通键盘减号渲染

# =========================================================================
# 1. 纯理想电路硬件参数（转换为 mpmath 的高精度浮点数 mpf）
# =========================================================================
N_cells = 20        # 晶胞段数
C1 = mpmath.mpf('1e-10')            # 胞内极弱耦合电容: 10 pF (高频贴片瓷片)
C2 = mpmath.mpf('1e-9')              # 胞间强耦合电容: 1 nF  (比C1大100倍，达成完美衰减)
L  = mpmath.mpf('8e-7')              # 节点接地并联电感: ~25.33 uH

Vs = mpmath.mpf('1.0')               # 交流信号源电动势幅值 (1 V)
Rs = mpmath.mpf('50.0')              # 工业标准射频源标准内阻 (50 欧姆)
R_L = mpmath.mpf('0')              # 电感电阻  (现在理论推导上完全不能允许有一点点电阻 有一点都不行 直接全部漏)

N_nodes = 2 * N_cells  # 全网总节点数（80个节点）

# =========================================================================
# 2. 计算纯电路的并联谐振频率 f0
# =========================================================================
# 在该频率下，各节点的对地电感导纳与电容总导纳(C1+C2)发生并联谐振，整体对地开路
omega0 = 1 / mpmath.sqrt(L * (C1 + C2))
f0 = omega0 / (2 * mpmath.mpf('2') * mpmath.pi)

# 纯电路理论手算分压校验值
z_in_theory = 1.0 / (omega0 * C2)                               # 输入等效电容阻抗: 10.493 欧姆
v0_theory = Vs * z_in_theory / mpmath.sqrt(Rs**2 + z_in_theory**2)    # 始端节点0的分压幅值: ~0.205 V

print(f"=== 纯电路稳态分析参数 ===")
print(f"节点并联谐振频率 f0 = {float(f0):.3f} Hz")
print(f"始端节点0理论分压值 |V0| = {float(v0_theory):.4f} V")

# =========================================================================
# 3. 构造基尔霍夫节点导纳矩阵 Y(ω)
# =========================================================================
def build_Y(omega):
    # 导纳就是阻抗导数
    # 创建 mpmath 专属的高精度 80x80 空矩阵
    Y = mpmath.matrix(N_nodes, N_nodes)

    YL  = 1 / (mpmath.j * omega * L + R_L)   # 电感导纳
    YC1 = mpmath.j * omega * C1        # 胞内电容导纳
    YC2 = mpmath.j * omega * C2        # 胞间电容导纳

    # 导纳矩阵的主对角线元素 Y[i, i]代表 所有直接连接在节点i上的支路导纳的总和
    # 各节点通用的 KCL 括号合并展开结构如下：
    # $$I_{\text{注入1}} = (Y_L + Y_{C1} + Y_{C2}) V_1 - Y_{C1} V_0 - Y_{C2} V_2$$

    # KCL项：各格点并联接地电感
    for i in range(N_nodes):
        Y[i, i] += YL

    # KCL项：跨接在 A, B 之间的胞内电容 C1
    for n in range(N_cells):
        A, B = 2*n, 2*n+1
        Y[A, A] += YC1; Y[B, B] += YC1
        # 这个负号是方程给出来的（源于电流取决于两端电压差，展开移项后邻居电压项前天然带减号）
        Y[A, B] -= YC1; Y[B, A] -= YC1

    # KCL项：跨接在当前 B 与下一级 A_next 之间的胞间电容 C2
    for n in range(N_cells - 1):
        B = 2*n + 1
        A_next = 2*(n+1)
        Y[B, B] += YC2; Y[A_next, A_next] += YC2
        # 这个负号同样由 KCL 邻居分流方程决定
        Y[B, A_next] -= YC2; Y[A_next, B] -= YC2

    return Y

# =========================================================================
# 4. 求解谐振点 f0 处的全网受迫交流稳态响应
# =========================================================================
Y = build_Y(omega0)

# 并在始端节点加上戴维南 电源内阻分压
Y[0, 0] += 1 / Rs

# 给最后一个79号接地 就是给他加一个爆大的导纳
Y[N_nodes-1, N_nodes-1] += 0

# 将Norton电流源并联到节点0
I = mpmath.matrix(N_nodes, 1)
I[0] = Vs / Rs

# 求解全网交流对地电压复数向量
# solve 就是解方程 YV=I (采用 mpmath 高精度 LU 分解盲解)
V = mpmath.lu_solve(Y, I)

# 计算高精度幅值，并转回普通的 float 供 matplotlib 画图
Vmag = [float(mpmath.fabs(v)) for v in V]

# 后台打印第1号节点（第一个B节点）的盲解状态，用底层数据检验“纯代数消元逼出的虚地”
print("\n=== 矩阵盲解节点状态校验 ===")
print(f"第一个节点 (Node 0) 实测电压幅值: {Vmag[0]:.6f} V")
print(f"第二个节点 (Node 1) 真实复数电压: {V[1]}")
print(f"第二个节点 (Node 1) 实测电压幅值: {Vmag[1]}")

# =========================================================================
# 5. 可视化：超深探测对数轴（彻底消灭坐标轴乱码，让波谷完美现形）
# =========================================================================
plt.figure(figsize=(12, 6))
node_indices = np.arange(len(V))

# -------------------------------------------------------------------------
# 🛠️ 核心修正：在 mpmath 高精度下直接提取 log10，彻底规避 float 转换导致的乱码
# -------------------------------------------------------------------------
# 1. 直接在高精度下把每个节点的对数电位算出来（限制最低到 -350，防止 log10(0) 报错）
log_Vmag = []
for v in V:
    mag = mpmath.fabs(v)
    if mag < mpmath.mpf('1e-350'):
        log_Vmag.append(-350.0)
    else:
        log_Vmag.append(float(mpmath.log10(mag)))

log_v0_theory = float(mpmath.log10(v0_theory))

# 2. 换成普通的 plt.plot。因为纵坐标已经是我们算好的指数了，所以画面天然就是完美的对数图！
plt.plot(node_indices, log_Vmag, 'o-', color='tab:blue', lw=1.2, ms=3, label='全网格点电压响应')
plt.plot(0, log_v0_theory, 'ro', markersize=8, zorder=5, label=f'手算等效单回路分压校验')

# 3. 强行卡死纵坐标。既然你现在的信号跌到了 10^-300，我们把底线斩到 -320
plt.ylim(-400, 10)

# 4. 手动生成干净、绝对不乱码的坐标轴刻度标签
tick_positions = np.arange(-300, 1, 50)  # 每隔 50 个数量级画一个刻度
tick_labels = [f'$10^{{{i}}}$' for i in tick_positions]
plt.yticks(tick_positions, tick_labels)

# -------------------------------------------------------------------------
# 标签与网格配置
# -------------------------------------------------------------------------
plt.xlabel("全网格点连续索引 (Node Index)", fontsize=11)
plt.ylabel("对地电压幅值 |V| (V)", fontsize=11)
plt.title("200节点高精度 KCL 矩阵硬解 —— 底部虚地（波谷）完美重现图", fontsize=12)

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11, loc='upper right')

plt.tight_layout()
plt.show()