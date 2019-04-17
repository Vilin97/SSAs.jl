using DiffEqJump, DiffEqBiological, DiffEqJump

function run_tests(methods = ["RSSA", "RSSA_m", "DirectCR", "DirectCR_dict", "Direct", "NRM"])
    jump_prob = load_model("DNA", Direct())
    times = range(jump_prob.prob.tspan...,length=101)
    for method in methods
        run_sum = 0
        N=8000
        for i in 1:N run_sum += my_solve(jump_prob, method, times)[2][end][3] end
        avg = run_sum/N
        if abs(avg - 592.655) < 10
            println("$method passed test")
        else
            println("$method FAILED test")
        end
    end
end
