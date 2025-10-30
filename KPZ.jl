using ProgressMeter
using Base.Threads

include("SFS/src/jl/integrators.jl") 


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack u = con
    @unpack bplan, fplan, p2 = tools
    
    @. f.k = - im * sqrt(p2) * φ.k
    mul!(f.x, bplan, f.k)
    @. f.x = u * f.x^2 / 2
    mul!(f.k, fplan, f.x)
    nothing
end

function run(n, seed, folder)
    d       = 1
    T       = 1e-2

    N       = 2^10
    L       = N/2
    Δt      = 2e-3
    
    con     = (u = 50.,)

    N_save  = 2^10
    N_step  = N_save * 100
    equ     = 0

    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=ETD2)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f, :f3]
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    @showprogress for i in 1:equ ETD!(fields, tools, f!, con, i) end

    save_opt = (
        save_names=save_names, folder=folder, 
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false,
        SAVEDATA=(:tools,)
    )

    save_first(tools, con, fields, save_opt)
    @showprogress for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end


function local_run()
    num = 5
    n   = 96
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        save_info(numbers, folder)
        run(n, seed, folder)
    end
end

local_run()
