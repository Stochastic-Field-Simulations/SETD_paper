using ProgressMeter
using Base.Threads

include("SFS/src/jl/integrators.jl") 


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack r, u = con
    
    @. f.x = - (( r + u * φ.x^2 ) * φ.x)
    nothing 
end


function run(n, seed, folder)
    d       = 2
    T       = 1e-2
    N       = 2^8
    L       = N * 2
    TIME    = 5e6
    Δt      = 4e-1

    N_save  = 100
    N_step  = Int(TIME/Δt)
    @assert N_save<=N_step
    @assert N_step%N_save==0

    rc = -0.0167
    lnt = LinRange(-4, log(10, abs(rc)), n)[seed]
    r = 10^lnt + rc

    con     = (r = r, u = 1)    
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=ETD2)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f, :f3]
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    save_opt = (
        save_names=save_names, folder=folder, 
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=false, SAVECORR=true, N_write=1_000,
        SAVEDATA=(:tools,)
    )

    save_first(tools, con, fields, save_opt)
    for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end

function local_run()
    num = 6
    n   = 48
    
    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        run(n, seed, folder)
    end
end

local_run() # This takes some time...
