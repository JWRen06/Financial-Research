# 付费核验与经营资金：信息反馈何时减少启动资本

Paid Verification and Operating Capital: When Feedback Reduces Initial Capital

当前版本：2026-11-01。论文正文及证明见 [paper.pdf](paper.pdf)，可编辑源文件见 [source](source/main.tex)。

## 复现

Python 3.12。在独立环境中安装本篇依赖，然后从本篇目录运行：

```sh
python -m pip install -r reproduction/requirements.txt
python reproduction/run.py
```

一项有限报价策略及其输入证书重算。2026-11-01论文已将保证样本报价推广到任意固定正系数；数值材料仅核验原a_H=2、a_L=20实例，不能当作一般系数的数值验证。

复现输入和程序仅针对当前论文的核验任务。全部数据为合成或模型内确定性证书，不包含实盘账户或真实市场收益。代码运行通过不替代一般性证明、原创性核查或外部同行审阅。

目录只保留当前论文、可编辑源文件及必要复现内容，不收录旧稿、讲解稿、审稿过程或研究日志。编译LaTeX需要中文字体；源文件使用Windows字体设置，其他系统需调整字体。
