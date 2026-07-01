using JuMP, COPT

const MOI = JuMP.MOI

# 使用 direct_model 而不是 Model
model = direct_model(COPT.Optimizer())

# 添加变量和约束
@variable(model, x >= 0)
@variable(model, y >= 0)
@constraint(model, x + y >= 1)
@objective(model, Min, x + 2y)

# 现在 backend(model) 直接返回 COPT.Optimizer
optimizer = backend(model)
prob = optimizer.prob  # 这行现在可以正常工作了

# 设置tune参数
MOI.set(optimizer, MOI.RawOptimizerAttribute("TuneTimeLimit"), 60.0)

# 执行参数调优
COPT.COPT_Tune(prob)

# 获取调优结果数量
num_results = Ref{Cint}()
COPT.COPT_GetIntAttr(prob, "TuneResults", num_results)
println("调优结果数量: ", num_results[])

# 加载最优调优参数
COPT.COPT_LoadTuneParam(prob, 0)

# 使用调优后的参数进行优化
optimize!(model)

# 查看结果
println("目标值: ", objective_value(model))