using DiffEqBiological, DiffEqJump, Random, Parameters, RandomNumbers.Xorshifts

function cr_solve(jump_prob; rt0=nothing, dep_graph=nothing, times=nothing, use_dict = false)
    parameters = cr_do_setup(jump_prob, rt0, dep_graph, times, use_dict)
    return cr_run_simulation(parameters...)
end

function cr_do_setup(jump_prob, rt0=nothing, dep_graph=nothing, times=nothing, use_dict = false)
    @unpack massaction_jump = jump_prob
    @unpack u0 = jump_prob.prob
    N = length(u0)
    M = get_num_majumps(massaction_jump)
    t0 = jump_prob.prob.tspan[1]
    if rt0 === nothing
        if use_dict
            rt = empty_rt_dict()
        else
            rt = empty_rt_vec(M)
        end
    else rt = copy(rt0)
    end
    if dep_graph === nothing
        dep_graph = DiffEqJump.make_dependency_graph(N, massaction_jump) end
    if times === nothing
        times = [t0]
        states = [u0]
        save_all_jumps = true
    else
        states = [u0 for i in times]
        save_all_jumps = false
    end

    return jump_prob, rt, dep_graph, times, states, save_all_jumps
end

function cr_run_simulation(jump_prob, rt, dep_graph, times, states, save_all_jumps)
    state = copy(jump_prob.prob.u0) #initial state
    @unpack massaction_jump = jump_prob
    M = get_num_majumps(massaction_jump) #number of reactions
    N = length(state) #number of species
    propensities = [DiffEqJump.evalrxrate(state, k, massaction_jump) for k in 1:M]
    for k in 1:M
        my_insert!(rt,k,propensities[k]) end
    t0, t_end = jump_prob.prob.tspan #start and end times
    t = t0
    rx_total_rate = sum(values(rt.gsums))
    idx = [2]
    counter = 0
    s = Xoroshiro128Star()
    while t < t_end
        t += Random.randexp()/rx_total_rate
        l = DM_search(rx_total_rate,rt.gsums,rt.idvec,s)
        group = rt.groups[l]
        @unpack jidxs,numrates,gmax = group
        if l != neg_pow_of_2
            mu = r_sample(numrates,jidxs,propensities,gmax,s)
        else # DM on the 0-bin
            mu = DM_search(rt.gsums[l],[propensities[rx_id] for rx_id in jidxs], jidxs, s)
        end
        DiffEqJump.executerx!(state,mu,massaction_jump) # update the current state
        record_state!(states, save_all_jumps,state,times,t,idx)
        @unpack jtogroup,groups,gsums = rt
        for k in dep_graph[mu]
            old_propensity = propensities[k]
            propensities[k] = DiffEqJump.evalrxrate(state, k, massaction_jump)
            propensity = propensities[k]
            l = jtogroup[k][1]
            group = groups[l]
            @unpack gmin,gmax = group
            if gmin < propensity <= gmax
                gsums[l] += (propensity - old_propensity)
            else
                swap!(rt,k,old_propensity,propensity)
            end
            rx_total_rate += (propensity - old_propensity)
        end
    end
    fill_states!(states,save_all_jumps,state,times,idx)
    return times,states
end

function r_sample(numrates,jidxs,propensities,gmax,s)
    mu = 0
    while true
        r = rand(s)*numrates
        k = Int(ceil(r))
        mu = jidxs[k]
        r =  ceil(r) - r
        if r*gmax < propensities[mu]
            return mu
        end
    end
    return mu
end
