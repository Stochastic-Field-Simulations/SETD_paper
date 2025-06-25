using ProgressMeter
using Base.Threads

include("SFS/src/jl/integrators.jl") 


function f!(fields, con, tools)
    @unpack φ, f = fields
    @unpack r, u = con
    
    @. f.x = - (( r + u * φ.x^2 ) * φ.x)
    nothing 
end

function run1(n, seed, folder)
    d       = 1
    T       = 1e-1
    N       = 2^8
    L       = N/4
    TIME    = 1e4
    Δt      = 1e-2

    N_save  = 1000
    N_step  = Int(TIME/Δt)
    equ     = N_step
    @assert N_save<=N_step
    @assert N_step%N_save==0

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
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false
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

    N_save  = 1000
    N_step  = Int(TIME/Δt)
    equ     = N_step
    @assert N_save<=N_step
    @assert N_step%N_save==0

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
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false
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

    N_save  = 100
    N_step  = Int(TIME/Δt)
    equ     = N_step
    @assert N_save<=N_step
    @assert N_step%N_save==0

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
        N_step=N_step, N_save=N_save, t_start=now(),
        SAVEFIELD=true, SAVECORR=false, N_write=false
    )

    save_first(tools, con, fields, save_opt)
    @showprogress for i in 1:N_step
        ETD!(fields, tools, f!, con, i)
        check_and_save(fields, tools, i, save_opt)
    end
end

function local_run1()
    num = 1
    n   = 6
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD1/$num/$seed/"
        save_info(numbers, folder)
        run1(n, seed, folder)
    end
end

function local_run2()
    num = 2
    n   = 6
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD1/$num/$seed/"
        save_info(numbers, folder)
        run2(n, seed, folder)
    end
end


function local_run3()
    num = 4 # Obs
    n   = 2^10
    
    numbers = (n=n)

    @threads for seed in 1:n
        folder  = "data/SETD1/$num/$seed/"
        save_info(numbers, folder)
        run3(n, seed, folder)
    end
end


# local_run1()
# local_run2()
# local_run3()
