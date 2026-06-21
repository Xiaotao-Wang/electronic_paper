import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

m, n = 2, 1
phi = (1 + np.sqrt(5)) / 2

print("=" * 80)
print(f"(m,n)=({m},{n}) 原始网格求和（无任何近似/跳过）")
print(f"黄金分割 α = {phi:.10f}")
print("=" * 80)


def compute_integral_raw(omega2LC, N=500):
    """
    最原始的方式：
    均匀网格，逐点计算被积函数，直接求和
    不跳过任何点，不设阈值
    """
    k1 = np.linspace(-np.pi, np.pi, N)
    k2 = np.linspace(-np.pi, np.pi, N)
    dk1 = k1[1] - k1[0]
    dk2 = k2[1] - k2[0]

    integral_real = 0.0
    integral_imag = 0.0

    for i in range(N):
        for j in range(N):
            # 被积函数分母
            den = omega2LC * np.sin(k1[i] / 2) ** 2 - np.sin(k2[j] / 2) ** 2

            # 被积函数分子
            num_real = 1 - np.cos(m * k1[i] + n * k2[j])
            num_imag = -np.sin(m * k1[i] + n * k2[j])

            # 直接相除（即使是无穷大也保留）
            val_real = num_real / den
            val_imag = num_imag / den

            integral_real += val_real * dk1 * dk2
            integral_imag += val_imag * dk1 * dk2

    return integral_real, integral_imag


def compute_integral_raw_vectorized(omega2LC, N=500):
    """
    向量化版本（速度快，等价于上面的双重循环）
    """
    k1 = np.linspace(-np.pi, np.pi, N)
    k2 = np.linspace(-np.pi, np.pi, N)
    dk1 = k1[1] - k1[0]
    dk2 = k2[1] - k2[0]

    K1, K2 = np.meshgrid(k1, k2)

    den = omega2LC * np.sin(K1 / 2) ** 2 - np.sin(K2 / 2) ** 2
    num_real = 1 - np.cos(m * K1 + n * K2)
    num_imag = -np.sin(m * K1 + n * K2)

    # 直接除，不设阈值
    val_real = num_real / den
    val_imag = num_imag / den

    integral_real = np.sum(val_real) * dk1 * dk2
    integral_imag = np.sum(val_imag) * dk1 * dk2

    return integral_real, integral_imag


# ============================================
# 测试几个点
# ============================================
print("\n不同 ω²LC 下的积分（原始方法，N=500）:")
print("-" * 70)
print(f"{'ω²LC':<12} {'积分实部':<25} {'积分虚部':<25}")
print("-" * 70)

test_points = [0.5, 1.0, phi, 1.8, 2.0]
for omega2LC in test_points:
    int_real, int_imag = compute_integral_raw_vectorized(omega2LC, N=500)
    marker = " <<< α" if abs(omega2LC - phi) < 1e-10 else ""
    print(f"{omega2LC:<12.6f} {int_real:<25.15f} {int_imag:<25.15f}{marker}")

# ============================================
# 检查奇点附近的行为
# ============================================
print(f"\n奇点分析: ω²LC = α = {phi:.10f}")
print("-" * 70)

# 看看分母在哪些点接近0
N_check = 200
k1 = np.linspace(-np.pi, np.pi, N_check)
k2 = np.linspace(-np.pi, np.pi, N_check)
K1, K2 = np.meshgrid(k1, k2)

den = phi * np.sin(K1 / 2) ** 2 - np.sin(K2 / 2) ** 2

# 找到分母最小的点
den_flat = den.flatten()
k1_flat = K1.flatten()
k2_flat = K2.flatten()

# 排序找最小的10个
idx_smallest = np.argsort(np.abs(den_flat))[:10]

print(f"分母最接近0的10个点:")
print(f"{'k1':<15} {'k2':<15} {'分母值':<20} {'分子虚部':<20} {'被积函数虚部':<20}")
print("-" * 90)

for idx in idx_smallest:
    k1_val = k1_flat[idx]
    k2_val = k2_flat[idx]
    den_val = den_flat[idx]
    num_val = -np.sin(2 * k1_val + k2_val)
    integrand_val = num_val / den_val
    print(f"{k1_val:<15.10f} {k2_val:<15.10f} {den_val:<20.15f} {num_val:<20.15f} {integrand_val:<20.15f}")

# ============================================
# 对称性验证
# ============================================
print(f"\n对称性验证（随机取点检查 f(k1,k2) + f(-k1,-k2) = 0）:")
print("-" * 70)

np.random.seed(12345)
test_k1 = np.random.uniform(-np.pi, np.pi, 10)
test_k2 = np.random.uniform(-np.pi, np.pi, 10)

print(f"{'k1':<12} {'k2':<12} {'f(k1,k2)':<25} {'f(-k1,-k2)':<25} {'和':<25}")
print("-" * 99)

for i in range(10):
    den_plus = phi * np.sin(test_k1[i] / 2) ** 2 - np.sin(test_k2[i] / 2) ** 2
    f_plus = -np.sin(2 * test_k1[i] + test_k2[i]) / den_plus

    den_minus = phi * np.sin(-test_k1[i] / 2) ** 2 - np.sin(-test_k2[i] / 2) ** 2
    f_minus = -np.sin(-2 * test_k1[i] - test_k2[i]) / den_minus

    print(f"{test_k1[i]:<12.6f} {test_k2[i]:<12.6f} {f_plus:<25.15f} {f_minus:<25.15f} {f_plus + f_minus:<25.15f}")

# ============================================
# 网格精度收敛测试
# ============================================
print(f"\n网格精度收敛测试 (ω²LC = α):")
print("-" * 70)
print(f"{'N':<8} {'积分实部':<25} {'积分虚部':<25}")
print("-" * 60)

prev_imag = None
for N in [100, 200, 300, 400, 500]:
    int_real, int_imag = compute_integral_raw_vectorized(phi, N=N)
    diff = abs(int_imag - prev_imag) if prev_imag is not None else 0
    print(f"{N:<8} {int_real:<25.15f} {int_imag:<25.15f} (Δ={diff:.2e})")
    prev_imag = int_imag

# ============================================
# 在 [1.6, 1.62] 区间密集采样
# ============================================
print(f"\n在 [1.6, 1.62] 区间密集采样（N=400）:")
print("-" * 70)

omega_values = np.sort(np.unique(np.concatenate([
    np.linspace(1.60, phi - 0.002, 10),
    np.linspace(phi - 0.002, phi - 0.0001, 20),
    [phi],
    np.linspace(phi + 0.0001, phi + 0.002, 20),
    np.linspace(phi + 0.002, 1.62, 10)
])))

results_real = []
results_imag = []

for i, om in enumerate(omega_values):
    if i % 10 == 0:
        print(f"进度: {i + 1}/{len(omega_values)}")
    int_r, int_i = compute_integral_raw_vectorized(om, N=400)
    results_real.append(int_r)
    results_imag.append(int_i)

results_real = np.array(results_real)
results_imag = np.array(results_imag)

# ============================================
# 绘图
# ============================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f'(m,n)=({m},{n}) 原始网格求和 N=400  α = {phi:.6f}',
             fontsize=14, fontweight='bold')

# 图1: 虚部全区间
ax = axes[0, 0]
ax.plot(omega_values, results_imag, 'b-o', markersize=4, linewidth=1.5)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('积分虚部', fontsize=12)
ax.set_title(f'虚部 [{omega_values[0]:.3f}, {omega_values[-1]:.3f}]', fontsize=13)
ax.grid(True, alpha=0.3)

# 图2: 虚部放大 α 附近
ax = axes[0, 1]
mask = np.abs(omega_values - phi) < 0.003
ax.plot(omega_values[mask], results_imag[mask], 'b-o', markersize=6, linewidth=2)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('积分虚部', fontsize=12)
ax.set_title(f'α附近放大 [{phi - 0.003:.6f}, {phi + 0.003:.6f}]', fontsize=13)
ax.grid(True, alpha=0.3)

# 图3: 实部全区间
ax = axes[1, 0]
ax.plot(omega_values, results_real, 'purple', linewidth=1.5)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('积分实部', fontsize=12)
ax.set_title('实部', fontsize=13)
ax.grid(True, alpha=0.3)

# 图4: 虚部数值表格
ax = axes[1, 1]
ax.axis('off')

mask_table = np.abs(omega_values - phi) < 0.003
table_text = f"α = {phi:.10f}\n\n"
table_text += "ω²LC              虚部\n"
table_text += "-" * 42 + "\n"
for i in np.where(mask_table)[0]:
    delta = omega_values[i] - phi
    marker = " <<< α" if abs(delta) < 1e-10 else ""
    table_text += f"{omega_values[i]:.10f}  {results_imag[i]:.12f}{marker}\n"

ax.text(0.1, 0.95, table_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.show()

# 最终结论
idx_phi = np.argmin(np.abs(omega_values - phi))
print(f"\n{'=' * 80}")
print(f"最终结果")
print(f"{'=' * 80}")
print(f"\n在 ω²LC = α = {phi:.10f}:")
print(f"  原始求和 实部 = {results_real[idx_phi]:.15f}")
print(f"  原始求和 虚部 = {results_imag[idx_phi]:.15f}")
print(f"\n理论分析:")
print(f"  被积函数虚部是奇函数: f(-k1,-k2) = -f(k1,k2)")
print(f"  在对称区间 [-π,π]×[-π,π] 上积分恒为 0")
print(f"  任何非零结果都是数值误差")