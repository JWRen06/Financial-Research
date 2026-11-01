# 未知报告平滑下的交易规则识别

Identifying Trading Rules under Unknown Reporting Smoothing

当前版本：2026-10-31。论文正文及证明见 [paper.pdf](paper.pdf)，可编辑源文件见 [source](source/paper.md)。

## 复现

Python 3.12。在独立环境中安装本篇依赖，然后从本篇目录运行：

```sh
python -m pip install -r reproduction/requirements.txt
python reproduction/run_experiments.py
python reproduction/continuous_experiments.py
```

原70组检查与连续扫描96组配置；这些配置不是独立随机样本。直接运行会在复现目录写结果，建议在副本执行。数值依赖为NumPy 2.3.5。

复现输入和程序仅针对当前论文的核验任务。全部数据为合成或模型内确定性证书，不包含实盘账户或真实市场收益。代码运行通过不替代一般性证明、原创性核查或外部同行审阅。

目录只保留当前论文、可编辑源文件及必要复现内容，不收录旧稿、讲解稿、审稿过程或研究日志。编译LaTeX需要中文字体；源文件使用Windows字体设置，其他系统需调整字体。
