import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# データ生成
x = np.linspace(0, 10, 100)
y = np.linspace(5, 0, 50)
data = np.random.random((len(y), len(x)))

# ヒートマップの作成
plt.figure(figsize=(10, 6))
sns.heatmap(data, xticklabels=x, yticklabels=y, cmap="viridis", alpha=0.7)

# 関数を重ねる (例: sin 関数)
x_vals = np.linspace(0, 10, 500)  # プロット用のX軸
y_vals = 2.5 + np.sin(x_vals)     # プロット用のY軸 (例: y = 2.5 + sin(x))
plt.plot(x_vals, y_vals, color="red", linewidth=2, label="y = 2.5 + sin(x)")

# 軸の設定
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Heatmap with Overlaid Function")
plt.legend()
plt.show()
