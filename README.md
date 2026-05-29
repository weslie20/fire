# Fire

This is a collection of Python scripts and Jupyter notebooks for machine learning / deep learning experiments.

## Git Push 速查

```bash
# 首次推送（新建仓库）
git remote add origin <仓库地址>
git push -u origin main

# 日常推送
git add .
git commit -m "提交说明"
git push

# 强制推送（谨慎使用）
git push --force

# 推送指定分支
git push origin <分支名>

# 删除远程分支
git push origin --delete <分支名>
```

## 文件说明

| 文件 | 说明 |
|------|------|
| `fire.py` | 主入口脚本 |
| `1.py` ~ `4.py` | 练习/实验脚本 |
| `demo.py` | 演示脚本 |
| `boston_housing_tutorial.ipynb` | 波士顿房价预测教程 |
| `pinn_heat_equation.ipynb` | PINN 热方程教程 |
