module SSAs

export my_solve

include("utils.jl")
include("CR_Structures.jl")
include("RSSA.jl")
include("DirectCR.jl")
include("Direct.jl")
include("NRM.jl")
include("benchmark.jl")
include("test.jl")

function my_solve(jump_prob, method::String, times=nothing)
    if method == "RSSA"
        times, states = rssa_solve(jump_prob, times=times)
    elseif method == "RSSA_m"
        times, states = rssa_solve(jump_prob, times=times, try_new = true)
    elseif method == "DirectCR"
        times, states = cr_solve(jump_prob, times=times)
    elseif method == "DirectCR_dict"
        times, states = cr_solve(jump_prob, times=times, use_dict = true)
    elseif method == "Direct"
        times, states = dm_solve(jump_prob, times=times)
    elseif method == "NRM"
        times, states = nrm_solve(jump_prob, times=times)
    end
    try return times, states
    catch
        println("Wrong method!") end
end

end # module
