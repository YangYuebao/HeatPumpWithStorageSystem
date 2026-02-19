# HeatPumpWithStorageSystem

## 测试案例说明

设计优化的测试案例位于`./calculations/situation22.jl`。

命令行使用流程：

1. 安装好julia并配置环境变量。本程序使用1.9.4版本；
2. 下载本项目的；
3. 在项目路径下启动CMD，按照5个线程启动julia：`julia --threads 5`
4. 运行`using Pkg;Pkg.activate(".");Pkg.instantiate()`安装项目依赖
5. 主要参数说明见`./calculations/situation22.jl`中注释。
6. 运行`include(joinpath(pwd(),"calculations","situation22.jl"))`来执行计算脚本


