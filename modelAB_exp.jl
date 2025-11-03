using ProgressMeter
using Base.Threads

include("SFS/src/jl/integrators.jl") 


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack r, u = con
    
    @. f.x = - ( r + u * φ.x^2 ) * φ.x
    
    mul!(f.k, tools.fplan, f.x)
end

function run(n, seed, folder)
    d       = 2
    T       = 1e-2
    N       = 2^8
    L       = N*8
    TIME    = 1e7
    Δt      = .5
    N_step  = Int(TIME/Δt)
    rc      = -0.0103 # must be less than actual rc, due to finite size
    lnt     = LinRange(-4.5, log(10, abs(rc)), n)[seed]

    con     = (r = 10^lnt + rc, u = 1.)    
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=ETD2)
    
    field_names = [:φ, :ξ, :f, :f3]
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    save_opt    = (
        save_names=save_names, folder=folder, SAVEFIELD=false, SAVECORR=true,
        N_step=N_step, N_save=100, N_write=1_000, SAVEDATA=(:tools,), t_start=now(),
    )

    save_first(tools, con, fields, save_opt)
    for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end

function run()
    num = 5
    n   = 48
    
    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        run(n, seed, folder)
    end
end

run() # This takes some time...
