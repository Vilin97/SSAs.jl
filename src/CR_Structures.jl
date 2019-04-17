# dynamic table data structure to store and update rates

using Parameters
using Printf

const neg_pow_of_2 = -100 #negligible power of 2

abstract type RateTable end

# one group of rates
mutable struct RateGroup{I,F}
gid::I                  #ID of the group
gmin::F
gmax::F
numrates::Int           #number of reactions in the group
jidxs::Vector{I}#reaction IDs
end

# table of groups
mutable struct RateTableDict{F,G,I,T} <: RateTable #F = float, G = group, I = int, T = tuple
minrate::F
maxrate::F
groups::Dict{I,G} #group ID to group
gsums::Dict{I,F} #sums of propensities for each group
jtogroup::Dict{I,T} #reaction ID to group ID. Dict to swap reactions around between groups
idvec::Vector{I} #array of group IDs sorted in descending order
end

mutable struct RateTableVec{F,G,I,T} <: RateTable #F = float, G = group, I = int, T = tuple
minrate::F
maxrate::F
groups::Dict{I,G} #group ID to group
gsums::Dict{I,F} #sums of propensities for each group
jtogroup::Vector{T} #reaction ID to group ID. Dict to swap reactions around between groups
idvec::Vector{I} #array of group IDs sorted in descending order
end

function empty_group(l)
    return RateGroup(l,2^(float(l-1)),2^float(l),0, Vector{Int}(undef,10))
end

function empty_rt_dict()
    return RateTableDict(0.0,0.0,Dict{Int,RateGroup{Int,Float64}}(),Dict{Int,Float64}(),Dict{Int,Tuple{Int,Int}}(),Int[]) end

function empty_rt_vec(M)
    return RateTableVec(0.0,0.0,Dict{Int,RateGroup{Int,Float64}}(),Dict{Int,Float64}(),Vector{Tuple{Int,Int}}(undef,M),Int[]) end

function clear!(group::RateGroup)
    group.numrates = 0
end

function clear!(rt::RateTableDict)
    for group in values(rt.groups)
        clear!(group)
    end
    for id in keys(rt.gsums)
        rt.gsums[id] = 0.0
    end
    rt.jtogroup = Dict{Int,Tuple{Int,Int}}()
end

function clear!(rt::RateTableVec)
    for group in values(rt.groups)
        clear!(group)
    end
    for id in keys(rt.gsums)
        rt.gsums[id] = 0.0
    end
        rt.jtogroup = Vector{Tuple{Int,Int}}(undef,length(rt.jtogroup))
end


function my_insert!(group,jumpid)
    @unpack numrates,jidxs = group
    numrates += 1
    if numrates > length(jidxs)
        push!(jidxs,jumpid)
    else
        jidxs[numrates] = jumpid
    end
    group.numrates = numrates
end

function update_jtogroup!(rt::RateTableVec,jumpid,l,numrates)
    @unpack jtogroup = rt
    if jumpid <= length(jtogroup)
        jtogroup[jumpid] = (l,numrates)
    else push!(jtogroup, (l,numrates))
    end
end

function update_jtogroup!(rt::RateTableDict,jumpid,l,numrates)
    @unpack jtogroup = rt
    jtogroup[jumpid] = (l,numrates)
end

@inline function my_insert!(rt::RateTable, jumpid, rate)
    if rate > 0
        l = ceil(Int, log(rate)/log(2))
    else l = neg_pow_of_2
    end
    if l in keys(rt.groups)
        group = rt.groups[l]
        @unpack gsums = rt
        my_insert!(group,jumpid)
        @unpack numrates = group
        update_jtogroup!(rt,jumpid,l,numrates)
        gsums[l] += rate
    else
        group = empty_group(l)
        rt.minrate = min(rt.minrate,2^float(l-1))
        if rt.minrate == 0.0
            rt.minrate = 2^float(l-1)
        end
        rt.maxrate = max(rt.maxrate,2^float(l))
        rt.groups[l] = group
        rt.gsums[l] = 0.0
        push!(rt.idvec,l)
        sort!(rt.idvec, rev=true, alg=InsertionSort)
        my_insert!(rt,jumpid,rate)
    end
end

@inline function swap!(rt::RateTable, jumpid, rate, new_rate)
    @unpack jtogroup,gsums = rt
    l,ind = jtogroup[jumpid]
    group = rt.groups[l]
    @unpack numrates,jidxs = group
    id_of_last_rx = jidxs[numrates] #id of the reaction that gets swapped out
    jidxs[ind] = id_of_last_rx #put the last rx to where the deleted rx was
    jtogroup[id_of_last_rx] = (l,ind) #update jtogroup accordingly
    numrates -= 1
    gsums[l] -= rate
    group.numrates = numrates
    my_insert!(rt,jumpid,new_rate)
end

function Base.show(io::IO, group::RateGroup)
    println("ID: $(group.gid), group range: ($(@sprintf("%.2f", group.gmin)),$(@sprintf("%.2f", group.gmax))), numrxs: $(group.numrates), reactions in group: $(group.jidxs[1:group.numrates])")
end

function Base.show(io::IO, rt::RateTable)
    println("Range: ($(rt.minrate),$(rt.maxrate)) \nRates per group: $(values(rt.gsums)) \nRxs to groups: $(rt.jtogroup) \nGroup IDs: $(rt.idvec)\n$(values(rt.groups)) ")
end
