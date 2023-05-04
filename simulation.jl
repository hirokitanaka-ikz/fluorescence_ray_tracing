@everywhere include("MonteCarloModel.jl")


@everywhere function process(crystal, beam, Δd, max_ref_count, pump_depletion, n_cycle::Int64)::Tuple
    df = DataFrame(λgen=Float64[], λesc=Float64[], ref_N=Int64[], x0=Float64[], y0=Float64[], z0=Float64[],
    x=Float64[], y=Float64[], z=Float64[], kx=Float64[], ky=Float64[], kz=Float64[])

    i::Int64 = 0
    error_counter::Int64 = 0

    while i < n_cycle
        ray = generate_ray(crystal, beam, pump_depletion)   # generate a ray
        no_error = ray_tracing!(ray, crystal, Δd, max_ref_count)   # ray tracing until ray escapes (no_error == true if no error)
        if no_error
            push!(
                df, (ray.λgen, ray.λ, ray.Nr, ray.r0[1], ray.r0[2], ray.r0[3],
                ray.r[1], ray.r[2], ray.r[3], ray.k[1], ray.k[2], ray.k[3])
            )
            i += 1
        else
            error_counter += 1
        end
    end
    return df, error_counter
end


function simulation(params)
    crystal = create_uniaxial_crystal(params["crystal"])
    beam = create_beam(params["beam"], crystal)
    N = params["simulation"]["N"]
    Δd = params["simulation"]["step"]
    if params["simulation"]["max_ref_count"] == 0
        max_ref_count = Inf
    else
        max_ref_count = params["simulation"]["max_ref_count"]
    end
    pump_depletion::Bool = params["simulation"]["pump_depletion"]

    # parallel processing
    nproc = length(procs())
    if N % nproc == 0
        n_cycles = ones(Int64, nproc) * Int64(N / nproc)
    else
        rest::Int64 = N % nproc
        n_cycles = ones(Int64, nproc) * div(N, nproc)
        n_cycles[1] += rest
    end
    fs = []     # vector to store Future objects from processes
    for m in n_cycles
        f = @spawn process(crystal, beam, Δd, max_ref_count, pump_depletion, m)
        push!(fs, f)
    end
    println("*****Simulation start: $nproc process(es)*****")
    df_all = DataFrame()
    Nerr::Int64 = 0
    for (i, f) in enumerate(fs)
        df, err_count = fetch(f)
        if i == 1
            df_all = df
        else
            df_all = vcat(df_all, df)
        end
        Nerr += err_count
    end
    println("total ray tracing error count: $Nerr / $N")
    println("*****Simulation end*****")
    return df_all
end
