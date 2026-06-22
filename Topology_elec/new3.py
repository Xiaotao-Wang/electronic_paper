import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

m, n = 1, 1
phi = 1
L, C = 1.0, 1.0

print("=" * 80)
print(f"(m,n)=({m},{n}) 完整阻抗 Z̃(ω) 计算")
print("Z̃(ω) = (-iωL/4π²) × ∫∫ (1-exp(i·(mk1+nk2)))/(ω²LC·sin²(k1/2)-sin²(k2/2)) dk1dk2")
print(f"黄金分割 α = {phi:.10f}")
print("=" * 80)


def compute_full_impedance(omega2LC, N=500):
    """
    计算完整的 Z̃(ω)
    """
    omega = np.sqrt(omega2LC / (L * C))

    k1 = np.linspace(-np.pi, np.pi, N)
    k2 = np.linspace(-np.pi, np.pi, N)
    dk1 = k1[1] - k1[0]
    dk2 = k2[1] - k2[0]

    K1, K2 = np.meshgrid(k1, k2)

    # 分母
    den = omega2LC * np.sin(K1 / 2) ** 2 - np.sin(K2 / 2) ** 2

    # 分子: 1 - exp(i*(2k1 + k2))
    num = 1 - np.exp(1j * (m * K1 + n * K2))

    # 被积函数（复数）
    integrand = num / den

    # 数值积分
    integral = np.sum(integrand) * dk1 * dk2

    # 乘以系数: -i*omega*L/(4*pi^2)
    coefficient = -1j * omega * L / (4 * np.pi ** 2)
    Z_tilde = coefficient * integral

    return Z_tilde, integral, omega


# 测试几个点
print("\n不同 ω²LC 下的完整阻抗:")
print("-" * 85)
print(f"{'ω²LC':<12} {'积分实部':<18} {'积分虚部':<18} {'Z̃实部':<18} {'Z̃虚部':<18}")
print("-" * 85)

test_points = [0.5, 1.0, phi, 1.8, 2.0]
for om2LC in test_points:
    Z, integral, omega = compute_full_impedance(om2LC, N=600)
    marker = " <<< α" if abs(om2LC - phi) < 1e-10 else ""
    print(f"{om2LC:<12.6f} {integral.real:<18.12f} {integral.imag:<18.12f} "
          f"{Z.real:<18.12f} {Z.imag:<18.12f}{marker}")

# ============================================
# 在 [1.6, 1.62] 密集采样
# ============================================
print(f"\n在 [1.6, 1.62] 密集采样 (N=500):")
print("-" * 70)

# 生成采样点
omega_values = np.sort(np.unique(np.concatenate([
    np.linspace(1.60, phi - 0.003, 15),
    np.linspace(phi - 0.003, phi - 0.0001, 30),
    [phi],
    np.linspace(phi + 0.0001, phi + 0.003, 30),
    np.linspace(phi + 0.003, 1.62, 15)
])))

Z_real_vals = []
Z_imag_vals = []
integral_real_vals = []
integral_imag_vals = []

for i, om in enumerate(omega_values):
    if (i + 1) % 20 == 0:
        print(f"进度: {i + 1}/{len(omega_values)}")

    Z, integral, omega = compute_full_impedance(om, N=500)

    Z_real_vals.append(Z.real)
    Z_imag_vals.append(Z.imag)
    integral_real_vals.append(integral.real)
    integral_imag_vals.append(integral.imag)

Z_real_vals = np.array(Z_real_vals)
Z_imag_vals = np.array(Z_imag_vals)
integral_real_vals = np.array(integral_real_vals)
integral_imag_vals = np.array(integral_imag_vals)

# ============================================
# 绘图
# ============================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle('(m,n)=(2,1) 完整阻抗 Z̃(ω)  L=C=1  α=' + f'{phi:.6f}',
             fontsize=14, fontweight='bold')

idx_phi = np.argmin(np.abs(omega_values - phi))

# 图1: 积分实部（非零！）
ax = axes[0, 0]
ax.plot(omega_values, integral_real_vals, 'b-', linewidth=1.5, alpha=0.8)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.plot(phi, integral_real_vals[idx_phi], 'ro', markersize=8)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('积分实部', fontsize=12)
ax.set_title('∫∫ Re[integrand] dk1dk2 ≠ 0', fontsize=13)
ax.grid(True, alpha=0.3)

# 图2: 积分虚部（应该≈0）
ax = axes[0, 1]
ax.plot(omega_values, integral_imag_vals, 'r-', linewidth=1.5, alpha=0.8)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.plot(phi, integral_imag_vals[idx_phi], 'ro', markersize=8)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('积分虚部', fontsize=12)
ax.set_title('∫∫ Im[integrand] dk1dk2 ≈ 0 (奇函数)', fontsize=13)
ax.grid(True, alpha=0.3)

# 图3: Z̃虚部（由积分实部决定）
ax = axes[0, 2]
ax.plot(omega_values, Z_imag_vals, 'g-', linewidth=2, alpha=0.8)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.plot(phi, Z_imag_vals[idx_phi], 'ro', markersize=8)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('Im[Z̃]', fontsize=12)
ax.set_title('Im[Z̃] = -ωL/(4π²) × 积分实部', fontsize=13)
ax.grid(True, alpha=0.3)

# 图4: Z̃虚部放大 α 附近
ax = axes[1, 0]
mask = np.abs(omega_values - phi) < 0.005
ax.plot(omega_values[mask], Z_imag_vals[mask], 'g-o', markersize=5, linewidth=2)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.plot(phi, Z_imag_vals[idx_phi], 'ro', markersize=10)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('Im[Z̃]', fontsize=12)
ax.set_title('Im[Z̃] 放大 [' + f'{phi - 0.005:.6f}' + ', ' + f'{phi + 0.005:.6f}' + ']', fontsize=13)
ax.grid(True, alpha=0.3)

# 图5: Z̃实部
ax = axes[1, 1]
ax.plot(omega_values, Z_real_vals, 'purple', linewidth=1.5, alpha=0.8)
ax.axvline(x=phi, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.set_xlabel('ω²LC', fontsize=12)
ax.set_ylabel('Re[Z̃]', fontsize=12)
ax.set_title('Re[Z̃] = ωL/(4π²) × 积分虚部 ≈ 0', fontsize=13)
ax.grid(True, alpha=0.3)

# 图6: 数据表格
ax = axes[1, 2]
ax.axis('off')

mask_table = np.abs(omega_values - phi) < 0.005
table_text = "α = " + f"{phi:.10f}" + "\n\n"
table_text += "ω²LC              Im[Z̃]\n"
table_text += "-" * 38 + "\n"
for i in np.where(mask_table)[0]:
    delta = omega_values[i] - phi
    marker = " <<< α" if abs(delta) < 1e-10 else ""
    table_text += f"{omega_values[i]:.10f} {Z_imag_vals[i]:.12f}{marker}\n"

ax.text(0.1, 0.95, table_text, transform=ax.transAxes, fontsize=7,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.show()

# 最终结果
print(f"\n{'=' * 80}")
print(f"最终结果")
print(f"{'=' * 80}")
print(f"\n在 ω²LC = α = {phi:.10f}:")
print(f"  二重积分实部: {integral_real_vals[idx_phi]:.15f}  ← 核心！不为零")
print(f"  二重积分虚部: {integral_imag_vals[idx_phi]:.15f}  ← 应该≈0")
print(f"  ω = {np.sqrt(phi):.10f}")
print(f"  Z̃实部: {Z_real_vals[idx_phi]:.15f}")
print(f"  Z̃虚部: {Z_imag_vals[idx_phi]:.15f}  ← 你要的值！")
print(f"\n公式: Z̃虚部 = -ωL/(4π²) × 积分实部")
print(f"     = {-np.sqrt(phi) / (4 * np.pi ** 2):.10f} × {integral_real_vals[idx_phi]:.10f}")
print(f"     = {Z_imag_vals[idx_phi]:.15f}")