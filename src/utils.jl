using ReactionNetworkImporters, DiffEqProblemLibrary.JumpProblemLibrary, RandomNumbers.Xorshifts, Plots, DiffEqJump

function get_dep_graph(jump_prob, massaction_jump, N)
    try jump_prob.discrete_jump_aggregation.dep_gr
    catch
        DiffEqJump.make_dependency_graph(N, massaction_jump) end
end

function record_state!(states, save_all_jumps,state,times,t,idx)
    if save_all_jumps
        push!(states, copy(state))
        push!(times, t)
    else
        temp = idx[1]
        while temp <= length(times) && times[temp] < t
            states[temp] = copy(state)
            temp += 1
        end
        idx[1] = temp
    end
end

function fill_states!(states,save_all_jumps,state,times,idx)
    if save_all_jumps
        temp = idx[1]
        while temp < length(times)
            states[idx] = state
            temp += 1
        end
    end
end

function DM_search(total_rate, propensities, idvec, s)
    r = total_rate*rand(s)
    run_sum = 0.0
    l = 0
    for k in 1:length(idvec) #choose group
        l = idvec[k] #l is group ID
        run_sum += propensities[l]
        if run_sum > r
            return l
        end
    end
    return l
end

# functions to load specific models

function dna_jump_prob(method)
    JumpProblemLibrary.importjumpproblems()
    prob = prob_jump_dnarepressor
    return JumpProblem(prob.discrete_prob,method,prob.network, save_positions = (false,false))
end

function multi_jump_prob(method)
    JumpProblemLibrary.importjumpproblems()
    prob = prob_jump_multistate
    return JumpProblem(prob.discrete_prob,method,prob.network, save_positions = (false,false))
end

function big_jump_prob(method)
    tf = 1.
    fname = "network_data\\blank_file_equilibration.net"
    prnbng = loadrxnetwork(BNGNetwork(), "BNGRepressilator", fname)
    rn = prnbng.rn
    u0 = round.(Int,prnbng.u₀)
    addjumps!(rn,build_regular_jumps=false, minimal_jumps=true)
    prob = DiscreteProblem(rn,u0, (0.,tf), prnbng.p)
    return JumpProblem(prob, method, rn, save_positions = (false,false))
end

function load_model(model_name::String, method)
    if model_name == "DNA"
        jump_prob = dna_jump_prob(method)
    elseif model_name == "B_cell"
        jump_prob = big_jump_prob(method)
    elseif model_name == "multistate"
        jump_prob = multi_jump_prob(method)
    end
    return jump_prob
end

# functions to plot benchmarking results

function get_median_times(results, model)
    methods = sort!(collect(keys(results[model])))
    times = zeros(Float64, length(methods))
    for k in 1:length(methods)
        times[k] = BenchmarkTools.median(results[model][methods[k]]).time/10^9 #get seconds
    end
    for k in 1:length(methods)
        name = methods[k]
        if length(name) >= 11 && name[1:11] == "DiffEqJump."
            methods[k] = name[12:end]
        end
    end
    return times, methods
end

function plot_results(results)
    plots = []
    models= collect(keys(results))
    for k in 1:length(models)
        times, methods = get_median_times(results, models[k])
        push!(plots,
        plot(xticks = (1:10,methods), times,
        linetype=:bar,
        yaxis="miliseconds", xaxis = (tickfontrotation = 30.0),
        label = "",
        title = models[k]) )
    end
    return plot(plots...)
end
