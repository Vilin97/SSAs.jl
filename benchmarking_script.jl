using BenchmarkTools, TimerOutputs

const to = TimerOutput()
reset_timer!(to)

@timeit to "generate B cell suite" big_suite = SSAs.make_big_suite()
@timeit to "run the benchmarks" big_results = run(big_suite, verbose = true, seconds = 30)
BenchmarkTools.save("benchmark results B cell.json", big_results)

show(to)
