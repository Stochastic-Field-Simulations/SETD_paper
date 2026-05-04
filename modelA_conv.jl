using ProgressMeter
using Base.Threads
using LinearAlgebra: mul!
using UnPack
using Dates

using SFS


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack r, u = con
    @unpack fplan = tools
    
    @. f.x = - (( r + u * φ.x^2 ) * φ.x)
    mul!(f.k, fplan, f.x)
end

function run_sim(Δt, m, seed, folder, time_step_type)
    d       = 1
    T       = 1e-4
    N       = 2^8
    L       = N/4
    TIME    = 1
    N_step  = Int(TIME/Δt)
    N_save  = 1

    con     = (r = .3, u = 3.,  bφ = 1.)
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=time_step_type)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f, :f3]
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    
    if seed==1; SAVEDATA=(:tools,)
    else; SAVEDATA = (:tools,); nothing end

    save_opt    = (
        save_names=save_names, folder=folder,
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false, SAVEDATA=SAVEDATA
    )
    SFS.init!(fields, tools, con)
    save_first_para_opt(tools, con, fields, save_opt; para=seed==1)
    for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end

function local_run()
    n   = 2^16 # ensemble size
    M   = 8
    Δt0 = 1/16

    nums    = (9, 10, 11)
    steps   = (ETD1, ETD2, IF)
    for m in 1:M
        Δt = Δt0 / (2^(m-1))
        println("Running m=$m, Δt=$Δt")
        @time @threads for seed in 1:n
            for i in eachindex(nums)
                folder  = "data/SETD_paper/$(nums[i])/$m/$seed/"
                run_sim(Δt, m, seed, folder, steps[i])
            end
        end
    end
end

local_run()
