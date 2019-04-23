using DiffEqBiological, DiffEqJump, Random, Parameters, RandomNumbers.Xorshifts

const spec_bound = 1.1
const r_spec_bound = 1/1.1

function rssa_solve(jump_prob; times=nothing, try_new = false)
    parameters = do_setup(jump_prob, times, try_new)
    return run_simulation(parameters...)
end


function do_setup(jump_prob, times, try_new)
    @unpack massaction_jump = jump_prob
    @unpack u0 = jump_prob.prob
    t0 = jump_prob.prob.tspan[1]
    N = length(u0)
    M = get_num_majumps(massaction_jump)
    dep_graph = get_dep_graph(jump_prob, massaction_jump, N)
    if times === nothing
        times = [t0]
        states = [u0]
        save_all_jumps = true
    else
        states = [u0 for i in times]
        save_all_jumps = false
    end
    Lstate = ceil.(Int,u0.*r_spec_bound)
    Hstate = ceil.(Int,u0.*spec_bound)
    Lprop = [DiffEqJump.evalrxrate(Lstate, k, massaction_jump) for k in 1:M]
    Hprop = [DiffEqJump.evalrxrate(Hstate, k, massaction_jump) for k in 1:M]
    return jump_prob, dep_graph, times, states, Lstate, Hstate, Lprop, Hprop, save_all_jumps, try_new
end

function run_simulation(jump_prob, dep_graph, times, states, Lstate, Hstate, Lprop, Hprop, save_all_jumps, try_new)
    state = copy(jump_prob.prob.u0)
    @unpack massaction_jump = jump_prob
    M = get_num_majumps(massaction_jump)
    N = length(state)
    t0, t_end = jump_prob.prob.tspan
    t = t0
    rx_total_rate = sum(Hprop)
    spec_to_dep_rxs = DiffEqJump.spec_to_dep_rxs_map(N, massaction_jump)
    idx = [2]
    if try_new
        props = [DiffEqJump.evalrxrate(state, k, massaction_jump) for k in 1:M]
        known = [true for k in 1:M]
    end
    s = Xoroshiro128Star()
    while t < t_end
        if try_new mu, t_inc = r_sample_new(rx_total_rate, Hprop, Lprop, M, state, massaction_jump, known, props,s)
        else mu, t_inc = r_sample(rx_total_rate, Hprop, Lprop, M, state, massaction_jump,s) end
        t += t_inc
        DiffEqJump.executerx!(state,mu,massaction_jump)
        record_state!(states, save_all_jumps,state,times,t,idx)
        for (spec_id,_) in massaction_jump.net_stoch[mu]
            ss = state[spec_id]
            if try_new
                for k in spec_to_dep_rxs[spec_id]
                     known[k] = false
                 end
            end
            if ss < Lstate[spec_id] || ss > Hstate[spec_id] #check if species is out of bound
                Lstate[spec_id] = ceil(Int,ss*r_spec_bound)
                Hstate[spec_id] = ceil(Int,ss*spec_bound)
                for k in spec_to_dep_rxs[spec_id]
                    Lprop[k] = DiffEqJump.evalrxrate(Lstate, k, massaction_jump)
                    new_prop_bound = DiffEqJump.evalrxrate(Hstate, k, massaction_jump)
                    rx_total_rate += (new_prop_bound - Hprop[k])
                    Hprop[k] = new_prop_bound
                end
            end
        end
    end
    fill_states!(states,save_all_jumps,state,times,idx)
    return times,states
end

function r_sample(rx_total_rate, Hprop, Lprop, M, state, massaction_jump, s)
    accepted = false
    mu = 0
    prod = 1.
    while !accepted
        mu = DM_search(rx_total_rate, Hprop, 1:M, s)
        r = rand(s)
        prod *= r
        temp = Hprop[mu]
        if r*temp <= Lprop[mu]
            accepted = true
        elseif r*temp <= DiffEqJump.evalrxrate(state, mu, massaction_jump)
            accepted = true
        end
    end
    return mu, -log(prod)/rx_total_rate
end

function r_sample_new(rx_total_rate, Hprop, Lprop, M, state, massaction_jump, known, props, s)
    accepted = false
    mu = 0
    prod = 1.
    while true
        mu = DM_search(rx_total_rate, Hprop, 1:M, s)
        r = rand(s)
        prod *= r
        temp = Hprop[mu]
        if r*temp <= Lprop[mu]
            return mu, -log(prod)/rx_total_rate
        elseif !known[mu]
            prop = DiffEqJump.evalrxrate(state, mu, massaction_jump)
            props[mu] = prop
            if r*temp <= prop return mu, -log(prod)/rx_total_rate end
        end
    end
    return mu, -log(prod)/rx_total_rate
end
