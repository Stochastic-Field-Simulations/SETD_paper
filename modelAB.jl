using ProgressMeter
using Base.Threads

include("SFS/src/jl/integrators.jl") 


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack r, u = con
    @unpack fplan = tools
    
    @. f.x = - (( r + u * φ.x^2 ) * φ.x)
    mul!(f.k, fplan, f.x)
end

function run1(n, seed, folder)
    d       = 1
    T       = 1e-1
    N       = 2^8
    L       = N/4
    TIME    = 1e5
    Δt      = 1e-2
    N_step  = Int(TIME/Δt)
    equ     = N_step

    r = 10^LinRange(.5, -3, n)[seed]
    con     = (r = r, u = 0.)
    
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=ETD1)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f] 
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    @showprogress for i in 1:equ ETD!(fields, tools, f!, con, i) end

    save_opt = (
        save_names=save_names, folder=folder, 
        N_step=N_step, N_save=1000, t_start=now(),
        SAVEFIELD=true, SAVECORR=true, N_write=false,
        SAVEDATA=(:tools,)
    )

    save_first(tools, con, fields, save_opt)
    @showprogress for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end


function run2(n, seed, folder)
    d       = 1
    T       = 1e-1
    N       = 2^8
    L       = N/4
    TIME    = 1e6
    Δt      = 1e-2
    N_step  = Int(TIME/Δt)
    equ     = N_step

    r = 10^LinRange(.5, -3, n)[seed]
    con     = (r = r, u = 0.)
    
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=true, time_step=ETD1)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f] 
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    @showprogress for i in 1:equ ETD!(fields, tools, f!, con, i) end

    save_opt = (
        save_names=save_names, folder=folder, 
        N_step=N_step, N_save=1000, t_start=now(),
        SAVEFIELD=true, SAVECORR=true, N_write=false,
        SAVEDATA=(:tools,)
    )

    save_first(tools, con, fields, save_opt)
    @showprogress for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end


function run3(n, seed, folder)
    d       = 1
    T       = 1e-1
    N       = 2^8
    L       = N
    TIME    = 2e1
    Δt      = 1e-2
    N_step  = Int(TIME/Δt)
    equ     = N_step

    con     = (r = 1., u = 1)
    
    sys     = System(d, N, L, Δt; T=T)
    tools   = Tools(sys; seed=seed, conserved=false, time_step=ETD2)

    # field, noise and aux fields
    field_names = [:φ, :ξ, :f, :f3]
    save_names  = [:φ,]
    fields      = get_fields(tools, field_names)

    @showprogress for i in 1:equ ETD!(fields, tools, f!, con, i) end

    save_opt = (
        save_names=save_names, folder=folder, 
        N_step=N_step, N_save=100, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false,
        SAVEDATA=(:tools,)
    )

    save_first(tools, con, fields, save_opt)
    @showprogress for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end

function run1()
    num = 1
    n   = 6
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        save_info(numbers, folder)
        run1(n, seed, folder)
    end
end

function run2()
    num = 2
    n   = 6
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        save_info(numbers, folder)
        run2(n, seed, folder)
    end
end


function run3()
    num = 3
    n   = 2^10
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD_paper/$num/$seed/"
        save_info(numbers, folder)
        run3(n, seed, folder)
    end
end


run1()
run2()
run3()
