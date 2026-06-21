import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib import rcParams
import warnings

warnings.filterwarnings('ignore')

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

# 黄金分割常数
phi = (1 + np.sqrt(5)) / 2  # ≈ 1.618
alpha = phi  # 黄金分割的根 α² + α - 1 = 0


def integrand_imag(k1, k2, m, n, omega2LC):
    """被积函数的虚部"""
    numerator = -np.sin(m * k1 + n * k2)
    denominator = omega2LC * np.sin(k1 / 2) ** 2 - np.sin(k2 / 2) ** 2

    # 处理分母接近零的情况
    mask = np.abs(denominator) < 1e-12
    denominator[mask] = np.sign(denominator[mask]) * 1e-12

    return numerator / denominator


def compute_impedance_imag(omega2LC, m=2, n=1, L=1.0, C=1.0, N=500):
    """
    计算阻抗的虚部
    """
    omega = np.sqrt(omega2LC / (L * C))

    # 在布里渊区上离散化
    k1 = np.linspace(-np.pi, np.pi, N)
    k2 = np.linspace(-np.pi, np.pi, N)
    dk = k1[1] - k1[0]

    K1, K2 = np.meshgrid(k1, k2)

    # 被积函数的虚部
    numerator_imag = -np.sin(m * K1 + n * K2)
    denominator = omega2LC * np.sin(K1 / 2) ** 2 - np.sin(K2 / 2) ** 2

    # 处理奇异点
    mask = np.abs(denominator) > 1e-10
    integrand = np.zeros_like(denominator)
    integrand[mask] = numerator_imag[mask] / denominator[mask]

    # 二维梯形积分
    integral_value = np.sum(integrand) * dk ** 2

    # 阻抗虚部
    Z_imag = omega * L / (4 * np.pi ** 2) * integral_value

    return Z_imag, omega


# 生成采样点
print("=" * 70)
print(f"参数设置: (m,n) = (2,1)，黄金分割 α² + α - 1 = 0，α = {phi:.10f}")
print("=" * 70)

# 关注点：α²LC 应该接近黄金分割的某个特征值
# 由于 (2,1) 对应黄金分割，我们可能关注 ω²LC = α ≈ 1.618 附近
# 同时也在 ω²LC = 1 附近密集采样（通用奇点）

golden_ratio = phi
print(f"\n黄金分割点: ω²LC = α = {golden_ratio:.6f}")
print(f"可能的关键点: ω²LC = 1 (分母零点), ω²LC = α ≈ {golden_ratio:.6f}")

# 生成密集采样点
# 区间1: 0到1（从下方趋近1）
below_1_uniform = np.linspace(0.01, 0.9, 15)
below_1_dense = 1 - np.logspace(-4, -0.05, 25)
below_1 = np.sort(np.unique(np.concatenate([below_1_uniform, below_1_dense])))

# 区间2: 1到2（从上方趋近1，并经过黄金分割点）
above_1_dense = 1 + np.logspace(-4, -0.05, 25)
above_1_mid = np.linspace(1.1, 2.0, 20)
# 在黄金分割点附近额外加密
golden_near = np.linspace(golden_ratio - 0.05, golden_ratio + 0.05, 20)
above_1 = np.sort(np.unique(np.concatenate([above_1_dense, above_1_mid, golden_near])))

# 合并所有点
omega2LC_all = np.sort(np.unique(np.concatenate([below_1, [1.0], above_1])))

print(f"\n采样策略:")
print(f"  区间 [0,1] 采样点数: {len(below_1)}")
print(f"  区间 [1,2] 采样点数: {len(above_1)}")
print(f"  总采样点数: {len(omega2LC_all)}")
print(f"\n关键点附近的采样:")
print(f"  在 ω²LC = 1 附近最小步长: {np.min(np.diff(omega2LC_all[np.abs(omega2LC_all - 1) < 0.1])):.2e}")
print(f"  在 ω²LC = α ≈ {golden_ratio:.4f} 附近的点数: {len(golden_near)}")

# 计算阻抗虚部
print(f"\n开始计算虚部 Im[Z̃(ω)]...")
m, n = 2, 1  # 黄金分割参数
L, C = 1.0, 1.0

results = {'omega2LC': [], 'Z_imag': []}

for i, omega2LC in enumerate(omega2LC_all):
    if (i + 1) % 15 == 0:
        print(f"进度: {i + 1}/{len(omega2LC_all)}")

    try:
        Z_imag, omega = compute_impedance_imag(omega2LC, m, n, L, C, N=500)

        results['omega2LC'].append(omega2LC)
        results['Z_imag'].append(Z_imag)

    except Exception as e:
        print(f"ω²LC = {omega2LC:.6f} 时计算失败: {e}")

results['omega2LC'] = np.array(results['omega2LC'])
results['Z_imag'] = np.array(results['Z_imag'])

print(f"成功计算 {len(results['omega2LC'])} 个点")

# 创建详细的分析图形
fig = plt.figure(figsize=(18, 12))
fig.suptitle(f'Im[Z̃(ω)] vs ω²LC  (m,n)=({m},{n})  黄金分割 α={phi:.4f}, α²+α-1=0',
             fontsize=14, fontweight='bold')

# 图1: 全区间概览 [0, 2] 标注黄金分割点
ax1 = plt.subplot(2, 3, 1)
ax1.plot(results['omega2LC'], results['Z_imag'], 'b-', linewidth=1.5, alpha=0.7)
ax1.axvline(x=1, color='red', linestyle='--', alpha=0.7, linewidth=1.5, label='ω²LC = 1')
ax1.axvline(x=golden_ratio, color='gold', linestyle='--', alpha=0.7, linewidth=1.5,
            label=f'ω²LC = α = {golden_ratio:.3f}')
ax1.set_xlabel('ω²LC', fontsize=12)
ax1.set_ylabel('Im[Z̃]', fontsize=12)
ax1.set_title('全区间概览 [0, 2]', fontsize=13)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax1.set_xlim(0, 2)

# 图2: 下方趋近1 [0.9, 1.0]
ax2 = plt.subplot(2, 3, 2)
mask_below = (results['omega2LC'] >= 0.9) & (results['omega2LC'] <= 1.0)
ax2.plot(results['omega2LC'][mask_below], results['Z_imag'][mask_below],
         'b-o', markersize=3, linewidth=1.5)
ax2.axvline(x=1, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
ax2.set_xlabel('ω²LC', fontsize=12)
ax2.set_ylabel('Im[Z̃]', fontsize=12)
ax2.set_title('从下方趋近1 [0.9, 1.0]', fontsize=13)
ax2.grid(True, alpha=0.3)

# 图3: 上方趋近1 [1.0, 1.1]
ax3 = plt.subplot(2, 3, 3)
mask_above = (results['omega2LC'] >= 1.0) & (results['omega2LC'] <= 1.1)
ax3.plot(results['omega2LC'][mask_above], results['Z_imag'][mask_above],
         'r-s', markersize=3, linewidth=1.5)
ax3.axvline(x=1, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
ax3.set_xlabel('ω²LC', fontsize=12)
ax3.set_ylabel('Im[Z̃]', fontsize=12)
ax3.set_title('从上方趋近1 [1.0, 1.1]', fontsize=13)
ax3.grid(True, alpha=0.3)

# 图4: 1附近放大 [0.99, 1.01]
ax4 = plt.subplot(2, 3, 4)
mask_near = (results['omega2LC'] >= 0.99) & (results['omega2LC'] <= 1.01)
mask_left = mask_near & (results['omega2LC'] < 1)
mask_right = mask_near & (results['omega2LC'] > 1)

if np.any(mask_left):
    ax4.plot(results['omega2LC'][mask_left], results['Z_imag'][mask_left],
             'b-o', markersize=4, linewidth=1.5, label='从下方')
if np.any(mask_right):
    ax4.plot(results['omega2LC'][mask_right], results['Z_imag'][mask_right],
             'r-s', markersize=4, linewidth=1.5, label='从上方')

ax4.axvline(x=1, color='green', linestyle='--', alpha=0.7, linewidth=1.5)
ax4.set_xlabel('ω²LC', fontsize=12)
ax4.set_ylabel('Im[Z̃]', fontsize=12)
ax4.set_title('1附近放大 [0.99, 1.01]', fontsize=13)
ax4.grid(True, alpha=0.3)
ax4.legend()

# 图5: 黄金分割点附近
ax5 = plt.subplot(2, 3, 5)
mask_golden = (results['omega2LC'] >= golden_ratio - 0.05) & (results['omega2LC'] <= golden_ratio + 0.05)
ax5.plot(results['omega2LC'][mask_golden], results['Z_imag'][mask_golden],
         'g-D', markersize=4, linewidth=1.5)
ax5.axvline(x=golden_ratio, color='gold', linestyle='--', alpha=0.7, linewidth=1.5,
            label=f'α = {golden_ratio:.4f}')
ax5.set_xlabel('ω²LC', fontsize=12)
ax5.set_ylabel('Im[Z̃]', fontsize=12)
ax5.set_title(f'黄金分割点 α 附近 [{golden_ratio - 0.05:.2f}, {golden_ratio + 0.05:.2f}]', fontsize=12)
ax5.grid(True, alpha=0.3)
ax5.legend()

# 图6: 虚部在三个关键区域的对比
ax6 = plt.subplot(2, 3, 6)
# 在0附近
mask_small = results['omega2LC'] < 0.1
# 在1附近
mask_near1 = np.abs(results['omega2LC'] - 1) < 0.05
# 在黄金分割附近
mask_near_phi = np.abs(results['omega2LC'] - golden_ratio) < 0.05

if np.any(mask_small):
    ax6.plot(results['omega2LC'][mask_small], results['Z_imag'][mask_small],
             'b.', markersize=6, label='ω²LC → 0', alpha=0.7)
if np.any(mask_near1):
    ax6.plot(results['omega2LC'][mask_near1], results['Z_imag'][mask_near1],
             'r.', markersize=6, label='ω²LC → 1', alpha=0.7)
if np.any(mask_near_phi):
    ax6.plot(results['omega2LC'][mask_near_phi], results['Z_imag'][mask_near_phi],
             'g.', markersize=6, label=f'ω²LC → α={golden_ratio:.3f}', alpha=0.7)

ax6.set_xlabel('ω²LC', fontsize=12)
ax6.set_ylabel('Im[Z̃]', fontsize=12)
ax6.set_title('三个关键区域对比', fontsize=13)
ax6.grid(True, alpha=0.3)
ax6.legend()

plt.tight_layout()
plt.show()

# 详细数据分析
print("\n" + "=" * 70)
print("详细数据分析")
print("=" * 70)

# ω²LC → 0 的行为
print(f"\n1. ω²LC → 0 时的行为 (前5个点):")
print(f"   {'ω²LC':<15} {'Im[Z̃]':<20}")
print(f"   " + "-" * 35)
for i in range(min(5, len(results['omega2LC']))):
    print(f"   {results['omega2LC'][i]:<15.10f} {results['Z_imag'][i]:<20.10f}")

# ω²LC → 1 的行为
print(f"\n2. ω²LC → 1 时的行为:")
below_1_final = results['omega2LC'] < 1
above_1_initial = results['omega2LC'] > 1

if np.any(below_1_final):
    print(f"\n   从下方趋近 (最后3个点):")
    indices = np.where(below_1_final)[0][-3:]
    for i in indices:
        delta = 1 - results['omega2LC'][i]
        print(f"   ω²LC = {results['omega2LC'][i]:.10f}, Im[Z̃] = {results['Z_imag'][i]:.10f}, Δ = {delta:.2e}")

if np.any(above_1_initial):
    print(f"\n   从上方趋近 (前3个点):")
    indices = np.where(above_1_initial)[0][:3]
    for i in indices:
        delta = results['omega2LC'][i] - 1
        print(f"   ω²LC = {results['omega2LC'][i]:.10f}, Im[Z̃] = {results['Z_imag'][i]:.10f}, Δ = {delta:.2e}")

# ω²LC → α (黄金分割) 的行为
print(f"\n3. ω²LC → α = {golden_ratio:.6f} 时的行为:")
mask_golden_near = np.abs(results['omega2LC'] - golden_ratio) < 0.01
if np.any(mask_golden_near):
    golden_points = np.where(mask_golden_near)[0]
    print(f"   附近的点:")
    for i in golden_points:
        delta = results['omega2LC'][i] - golden_ratio
        print(f"   ω²LC = {results['omega2LC'][i]:.10f}, Im[Z̃] = {results['Z_imag'][i]:.10f}, Δ = {delta:.2e}")

# 极限分析
print(f"\n4. 极限外推:")

if np.any(below_1_final):
    last_below_idx = np.where(below_1_final)[0][-3:]
    if len(last_below_idx) >= 2:
        x_below = results['omega2LC'][last_below_idx]
        y_below = results['Z_imag'][last_below_idx]
        coeffs = np.polyfit(x_below, y_below, 1)
        limit_below = np.polyval(coeffs, 1.0)
        print(f"   从下方 → 1: Im[Z̃] ≈ {limit_below:.10f}")

if np.any(above_1_initial):
    first_above_idx = np.where(above_1_initial)[0][:3]
    if len(first_above_idx) >= 2:
        x_above = results['omega2LC'][first_above_idx]
        y_above = results['Z_imag'][first_above_idx]
        coeffs = np.polyfit(x_above, y_above, 1)
        limit_above = np.polyval(coeffs, 1.0)
        print(f"   从上方 → 1: Im[Z̃] ≈ {limit_above:.10f}")

# 输出黄金分割验证
print(f"\n5. 黄金分割验证:")
print(f"   α = {phi:.10f}")
print(f"   α² = {phi ** 2:.10f}")
print(f"   α² + α - 1 = {phi ** 2 + phi - 1:.10f} (应该 ≈ 0)")