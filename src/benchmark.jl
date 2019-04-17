using BenchmarkTools, DiffEqJump, DiffEqBiological

function make_small_suite()
    methods = [RSSA(), NRM(), Direct(), DirectCR()]
    my_methods = ["RSSA", "RSSA_m", "DirectCR", "DirectCR_dict", "Direct", "NRM"]
    models = ["DNA", "multistate"]
    suite = BenchmarkGroup()
    for model_name in models
        suite[model_name] = BenchmarkGroup([])
        jump_prob = load_model(model_name, Direct())
        times = range(jump_prob.prob.tspan...,length=101)
        for method in my_methods
            suite[model_name][string("my ",method)] = @benchmarkable my_solve($jump_prob,$method,$times)
        end
        saveat = jump_prob.prob.tspan[2]/100.
        for method in methods
            jump_prob = load_model(model_name, method)
            suite[model_name][string(method)[1:end-2]] = @benchmarkable solve($jump_prob,$(SSAStepper()),saveat = $saveat)
        end
        println("Done with ", model_name)
    end
return suite
end

function make_big_suite()
    methods = [RSSA(), NRM(), DirectCR()]
    my_methods = ["RSSA", "RSSA_m", "NRM", "DirectCR", "DirectCR_dict"]
    models = ["B_cell"]
    suite = BenchmarkGroup()
    for model_name in models
        suite[model_name] = BenchmarkGroup([])
        jump_prob = load_model(model_name, Direct())
        times = range(jump_prob.prob.tspan...,length=101)
        for method in my_methods
            suite[model_name][string("my ",method)] = @benchmarkable my_solve($jump_prob,$method,$times) samples = 5 seconds = 30
        end
        saveat = jump_prob.prob.tspan[2]/100.
        for method in methods
            jump_prob = load_model(model_name, method)
            suite[model_name][string(method)[1:end-2]] = @benchmarkable solve($jump_prob,$(SSAStepper()),saveat = $saveat) samples = 5 seconds = 30
        end
        println("Done with ", model_name)
    end
return suite
end
