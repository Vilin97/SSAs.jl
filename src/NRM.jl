using DataStructures, Random, Parameters, RandomNumbers.Xorshifts,  DiffEqBiological, DiffEqJump

function nrm_solve(jump_prob; dep_graph=nothing, times=nothing)
    parameters = nrm_do_setup(jump_prob, dep_graph, times)
    return nrm_run_simulation(parameters...)
end


function nrm_do_setup(jump_prob, dep_graph, times)
    @unpack massaction_jump = jump_prob
    @unpack u0 = jump_prob.prob
    t0 = jump_prob.prob.tspan[1]
    N = length(u0)
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
    return jump_prob, dep_graph, times, states, save_all_jumps
end

function nrm_run_simulation(jump_prob, dep_graph, times, states, save_all_jumps)
    state = copy(jump_prob.prob.u0)
    @unpack massaction_jump = jump_prob
    M = get_num_majumps(massaction_jump)
    N = length(state)
    t0, t_end = jump_prob.prob.tspan
    t = t0
    idx = [2]
    propensities = [DiffEqJump.evalrxrate(state, k, massaction_jump) for k in 1:M]
    h = MutableBinaryMinHeap(Random.randexp(M)./propensities)
    s = Xoroshiro128Star()
    while t < t_end
        t, mu = top_with_handle(h)
        DiffEqJump.executerx!(state,mu,massaction_jump)
        record_state!(states, save_all_jumps, state, times, t, idx)
        for k in dep_graph[mu]
            old_propensity = propensities[k]
            propensities[k] = DiffEqJump.evalrxrate(state, k, massaction_jump)
            if old_propensity != 0
                update!(h, k, old_propensity/propensities[k]*(h[k]-t)+t)
            else
                update!(h, k, t+Random.randexp(s)/propensities[k])
            end
        end
        update!(h, mu, t+Random.randexp(s)/propensities[mu])
    end
    fill_states!(states,save_all_jumps,state,times,idx)
    return times, states
end
