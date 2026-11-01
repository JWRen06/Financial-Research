# 有限预算下的机会捕获：利润范围与测量精度

Opportunity Capture under Finite Budgets: Profit Bounds and Measurement Precision

当前版本：2026-10-31。论文正文及证明见 [paper.pdf](paper.pdf)，可编辑源文件见 [source](source/main.tex)。

## 复现

Python 3.12。在独立环境中安装本篇依赖，然后从本篇目录运行：

```sh
python -m pip install -r reproduction/requirements.txt
python reproduction/run_all.py
```

连续投入与有限奖励的合成案例。入口运行四个原程序、聚焦案例与停止规则核查，写入run_output；当前结果与更早环境的逐路径结果可能不同。

复现输入和程序仅针对当前论文的核验任务。全部数据为合成或模型内确定性证书，不包含实盘账户或真实市场收益。代码运行通过不替代一般性证明、原创性核查或外部同行审阅。

目录只保留当前论文、可编辑源文件及必要复现内容，不收录旧稿、讲解稿、审稿过程或研究日志。编译LaTeX需要中文字体；源文件使用Windows字体设置，其他系统需调整字体。
