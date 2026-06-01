# API 参考

## 核心模块

### transfer_matrix

薄膜 Transfer-Matrix Method 实现。

```python
from polarization_poincare import transfer_matrix

# 计算单层膜的 transfer matrix
M = transfer_matrix(n1, n2, d, wavelength, theta)
```

### mueller

Jones → Mueller 矩阵转换。

```python
from polarization_poincare.mueller import jones_to_mueller

M = jones_to_mueller(J)
# M: 4×4 Mueller matrix
```

### berreman

Berreman 4×4 求解器，用于各向异性膜层。

```python
from polarization_poincare.berreman import solve

# 求解 4×4 特征值问题
psi, delta = solve(epsilon_tensor, k, d)
```

### anisotropic

各向异性介质折射率张量处理。

```python
from polarization_poincare.anisotropic import RefractiveIndexTensor

tensor = RefractiveIndexTensor(n_o=1.5, n_e=1.6, axis=(0, 0, 1))
```

## 参数说明

| 参数 | 类型 | 说明 |
|------|------|------|
| `n1, n2` | complex | 入射/出射介质折射率 |
| `d` | float | 膜层厚度 (nm) |
| `wavelength` | float | 真空波长 (nm) |
| `theta` | float | 入射角 (deg) |