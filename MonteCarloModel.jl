using LinearAlgebra
using Distributions
using Interpolations
using CSV, DataFrames

include("Crystal.jl")
include("Beam.jl")
include("Ray.jl")
include("CosineCorrector.jl")

SMALL = 1e-8    # small value to avoid errors due to floating point


# load fluorescence and absorption spectra of an isotropic crystal
function load_spectra_isotropic(name::String, T::Int64)::Tuple
    filepath = "spectra/$(name)_$(T)K.csv"
    df = CSV.read(filepath, header=1, DataFrame)
    λ_vector::Vector{Float64} = df[:,1]
    If::Vector{Float64} = df[:,2]
    σabs::Vector{Float64} = df[:,3]
    return λ_vector, If, σabs
end


# load fluorescence and absorption spectra of an uniaxial crystal
function load_spectra_uniaxial(name::String, T::Int64)::Tuple
    filepath_π = "spectra/$(name)_pi_$(T)K.csv"
    filepath_σ = "spectra/$(name)_sigma_$(T)K.csv"
    df_π = CSV.read(filepath_π, header=1, DataFrame)
    df_σ = CSV.read(filepath_σ, header=1, DataFrame)
    λ_vector::Vector{Float64} = df_π[:,1]
    If_π::Vector{Float64} = df_π[:,2]
    If_σ::Vector{Float64} = df_σ[:,2]
    σabs_π::Vector{Float64} = df_π[:,3] # [cm^2]
    σabs_σ::Vector{Float64} = df_σ[:,3] # [cm^2]
    return λ_vector, If_π, If_σ, σabs_π, σabs_σ
end


function create_crystal(params_crystal)
    if params_crystal["num_axes"] == 1
        return create_isotropic_crystal(params_crystal)
    elseif params_crystal["num_axes"] == 2
        return create_uniaxial_crystal(params_crystal)
    # elseif params_crystal["num_axes"] == 3
    #     return create_biaxial_crystal(params_crystal)
    else
        prinln("Number of axes has to be 1, 2, or 3!")
    end
end


# define isotropic crystal
function create_isotropic_crystal(params_crystal)::Crystal_isotropic
    name = params_crystal["name"]
    T = params_crystal["T"]
    n = params_crystal["n1"]
    conc = params_crystal["cation_density"] * params_crystal["doping_level"] * params_crystal["correction_coeff"] / 100
    W = params_crystal["W"]
    H = params_crystal["H"]
    L = params_crystal["L"]
    QE = params_crystal["QE"]
    alpha_b = params_crystal["alpha_b"]
    if params_crystal["shape"] == "cuboid"
        θ = 0.5π
    elseif params_crystal["shape"] == "brewster"
        θ = atan(n)
    else
        println("crystal's shape is neither 'cuboid' nor 'brewster'!")
    end
    ϕ = params_crystal["angle"] / 180 * π
    λ_vector, If, σabs = load_spectra_isotropic(name, T)
    α::Vector{Float64} = σabs * conc    # [cm^-1]
    p_planes::Tuple = ([0, 0, 0], [W, 0, 0], [0, 0, 0], [0, H, 0], [0, 0, 0], [0, 0, L])
    plane_normals::Tuple = ([1, 0, 0], [-sin(ϕ), -cos(ϕ), 0], [0, 1, 0], [0, -1, 0], [-cos(θ), 0, sin(θ)], [cos(θ), 0, -sin(θ)])
    return Crystal_isotropic(W, H, L, n, θ, T, conc, λ_vector, If, α, p_planes, plane_normals, QE, alpha_b)
end

# define uniaxial crystal
function create_uniaxial_crystal(params_crystal)::Crystal_uniaxial
    name = params_crystal["name"]
    T = params_crystal["T"]
    no = params_crystal["n1"]
    ne = params_crystal["n2"]
    conc = params_crystal["cation_density"] * params_crystal["doping_level"] * params_crystal["correction_coeff"] / 100
    W = params_crystal["W"]
    H = params_crystal["H"]
    L = params_crystal["L"]
    caxis = params_crystal["caxis"]
    QE = params_crystal["QE"]
    alpha_b = params_crystal["alpha_b"]
    if params_crystal["shape"] == "cuboid"
        θ = 0.5π
    elseif params_crystal["shape"] == "brewster"
        θ = atan(ne)
    else
        println("crystal's shape is neither 'cuboid' nor 'brewster'!")
    end
    ϕ = params_crystal["angle"] / 180 * π
    λ_vector, If_π, If_σ, σabs_π, σabs_σ = load_spectra_uniaxial(name, T)
    α_π::Vector{Float64} = σabs_π * conc
    α_σ::Vector{Float64} = σabs_σ * conc
    p_planes::Tuple = ([0, 0, 0], [W, 0, 0], [0, 0, 0], [0, H, 0], [0, 0, 0], [0, 0, L])
    plane_normals::Tuple = ([1, 0, 0], [-sin(ϕ), -cos(ϕ), 0], [0, 1, 0], [0, -1, 0], [-cos(θ), 0, sin(θ)], [cos(θ), 0, -sin(θ)])
    return Crystal_uniaxial(W, H, L, caxis, no, ne, θ, T, conc, λ_vector, If_π, If_σ, α_π, α_σ, p_planes, plane_normals, QE, alpha_b)
end


# find intersections of a line with vector k passing through a point r with crystal's planes
function find_intersection(r::Vector{Float64}, k::Vector{Float64}, crystal)::Tuple
    t_list = [dot((crystal.p_planes[i] - r), crystal.plane_normals[i]) / dot(k, crystal.plane_normals[i]) for i in 1:6]
    intersections = [r + t * k for t in t_list if abs(t) != Inf]
    return intersections[1], intersections[2]
end


# define beam
function create_beam(params_beam, crystal)::Beam
    λ = params_beam["wl"]
    k = params_beam["k"]
    p0 = params_beam["p0"]
    E = params_beam["E"]
    p1, p2 = find_intersection(p0, k, crystal)
    path_length = norm(p2 - p1)
    d = params_beam["d"]
    distribution = params_beam["distribution"]
    return Beam(λ, normalize(k), d, p0, p1, p2, path_length, E, distribution)
end


# generate random vector (3D random direction)
function generate_random_vector()::Vector{Float64}
    v = randn(3)
    v /= norm(v)
    return v
end


# generate random vector perpendicular to n
function generate_perpendicular_vector(n::Vector{Float64})::Vector{Float64}
    m = randn(3)
    v = cross(n, m)
    v /= norm(v)
    return v
end


# generate random wavelength from spectrum
function get_fluorescence_wavelength(E::Vector{Float64}, crystal::Crystal_isotropic)::Float64
    λ = crystal.λ_vector
    spectrum = crystal.If
    spectrum /= maximum(spectrum)
    f = LinearInterpolation(λ, spectrum)
    while true
        λrand = rand(Uniform(λ[1], λ[end]))
        Ifrand = rand()
        if Ifrand < f(λrand)
            return λrand
        end
    end
end


# generate random wavelength from spectrum
function get_fluorescence_wavelength(E::Vector{Float64}, crystal::Crystal_uniaxial)::Float64
    θ = acos( abs( dot(E, crystal.caxis) ) )
    coeff_π = cos(θ) / (cos(θ) + sin(θ))
    coeff_σ = 1 - coeff_π
    λ = crystal.λ_vector
    spectrum = coeff_π * crystal.If_π + coeff_σ * crystal.If_σ
    spectrum /= maximum(spectrum)
    f = LinearInterpolation(λ, spectrum)
    while true
        λrand = rand(Uniform(λ[1], λ[end]))
        Ifrand = rand()
        if Ifrand < f(λrand)
            return λrand
        end
    end
end


function point_in_crystal(r::Vector{Float64}, crystal)::Bool
    for i in range(1, 6)
        p = crystal.p_planes[i]
        n = crystal.plane_normals[i]
        if dot([r[1] - p[1], r[2] - p[2], r[3] - p[3]], n) < SMALL
            return false
        end
    end
    return true
end


# function that rotates a vector which is normal to [0,0,1] to be perpendicular to an arbitrary unit vector k
function rotate_vector(r::Vector{Float64}, k::Vector{Float64})
    axis = cross([0, 0, 1], k)
    angle = acos(dot([0, 0, 1], k))
    if angle == 0
        return r
    end
    r_rotated = cos(angle) * r + sin(angle) * cross(axis, r) + (1 - cos(angle)) * dot(axis, r) * axis
    return r_rotated
end


function get_ray_position(crystal, beam::Beam, pump_depletion::Bool)::Vector{Float64}
    if pump_depletion
        # Beer-Lambert law
        α = get_absorption_coefficient(beam, crystal) * 0.1     # [mm^-1]
        val_rand = rand(Uniform(exp(-α * beam.path_length), 1.0))
        p = beam.p1 + beam.k * ( -log(val_rand) / α )
    else
        # uniform distribution
        p = beam.p1 + beam.k * beam.path_length * rand()
    end

    if beam.distribution == "gauss"
        R = abs(randn()) * beam.d / 2 * 0.5
    elseif beam.distribution == "tophat"
        R = sqrt(rand())
    else
        println("beam distribution is neither 'gauss' nor 'tophat'!")
    end

    θ = 2π * rand()
    r = p + [R * cos(θ), R * sin(θ), 0.0]
    r_rotated = rotate_vector(r, beam.k)
    return r_rotated
end


# generate ray
function generate_ray(crystal, beam::Beam, pump_depletion::Bool)::Ray
    k = generate_random_vector()
    E = generate_perpendicular_vector(k)
    λ = get_fluorescence_wavelength(E, crystal)
    ray_pos = get_ray_position(crystal, beam, pump_depletion)
    while !point_in_crystal(ray_pos, crystal)
        ray_pos = get_ray_position(crystal, beam, pump_depletion)
    end
    return Ray(ray_pos, k, E, λ)
end


# find intersect of ray with a plane
function find_intersection(ray::Ray, crystal)::Tuple
    t_list = [dot((crystal.p_planes[i] - ray.r), crystal.plane_normals[i]) / dot(ray.k, crystal.plane_normals[i]) for i in 1:6]
    t_list = [t > 1e-10 ? t : Inf for t in t_list]  # small value to avoid an error due to floating point
    tmin = minimum(t_list)
    normal = crystal.plane_normals[argmin(t_list)]
    return tmin, normal
end


# forward ray
function forward!(ray::Ray, Δd::Float64)
    ray.r += ray.k * Δd
end


# redirect ray
function redirect!(ray::Ray)
    ray.k = generate_random_vector()
    ray.E = generate_perpendicular_vector(ray.k)
end


# reset wavelength of ray
function reset_wavelength!(ray::Ray, crystal)
    ray.λ = get_fluorescence_wavelength(ray.E, crystal)
end


# reflection of ray at a plane perpendicular to vector normal
function reflect!(ray::Ray, normal::Vector{Float64})
    ray.k -= 2 * dot(normal, ray.k) * normal
    ray.Nr += 1     # number of reflection +1
end


# redirect the transmitted ray using Snell's law (isotropic crystal)
function transmit!(ray::Ray, crystal::Crystal_isotropic, normal::Vector{Float64})
    θi = acos(dot(-normal, ray.k) / (norm(normal) * norm(ray.k)))   # indicent angle
    if θi < SMALL
        return ki
    end
    θo = asin(crystal.n * sin(θi))          # angle of refraction (Snell's law)
    Δθ = θo - θi
    l = normalize( normal - dot(ray.k, normal) * ray.k )
    ray.k = ray.k * cos(Δθ) + l * sin(Δθ)
end


# redirect the transmitted ray using Snell's law (uniaxial crystal)
function transmit!(ray::Ray, crystal::Crystal_uniaxial, normal::Vector{Float64})
    θi = acos(dot(-normal, ray.k) / (norm(normal) * norm(ray.k)))   # indicent angle
    if θi < SMALL
        return ki
    end
    n = minimum([crystal.no, crystal.ne])
    θo = asin(n * sin(θi))          # angle of refraction (Snell's law) using the average refractive index
    Δθ = θo - θi
    l = normalize( normal - dot(ray.k, normal) * ray.k )
    ray.k = ray.k * cos(Δθ) + l * sin(Δθ)
end


# get absorption coefficient for fluorescence ray in cm^-1 (isotropic crystal)
function get_absorption_coefficient(ray::Ray, crystal::Crystal_isotropic)::Float64
    f = LinearInterpolation(crystal.λ_vector, crystal.α)
    return f(ray.λ)
end


# get absorption coefficient for fluorescence ray in cm^-1 (uniaxial crystal)
function get_absorption_coefficient(ray::Ray, crystal::Crystal_uniaxial)::Float64
    θ = acos( abs( dot(ray.E, crystal.caxis) ) )
    coeff_π = cos(θ) / (cos(θ) + sin(θ))
    coeff_σ = 1 - coeff_π
    α = coeff_π * crystal.α_π + coeff_σ * crystal.α_σ
    f = LinearInterpolation(crystal.λ_vector, α)
    return f(ray.λ)
end


# get absorption coefficient for excitation beam in cm^-1 (isotropic crystal)
function get_absorption_coefficient(beam, crystal::Crystal_isotropic)::Float64
    f = LinearInterpolation(crystal.λ_vector, crystal.α)
    return f(beam.λ)
end


# get absorption coefficient for excitation beam in cm^-1 (uniaxial crystal)
function get_absorption_coefficient(beam, crystal::Crystal_uniaxial)::Float64
    θ = acos( abs( dot(beam.E, crystal.caxis) ) )
    coeff_π = cos(θ) / (cos(θ) + sin(θ))
    coeff_σ = 1 - coeff_π
    α = coeff_π * crystal.α_π + coeff_σ * crystal.α_σ   # [cm^-1]
    f = LinearInterpolation(crystal.λ_vector, α)
    return f(beam.λ)
end


# judge if ray is absorbed (true if absorbed)
function judge_absorbed(ray::Ray, crystal, d::Float64)::Bool
    α = get_absorption_coefficient(ray, crystal)
    A = 1 - exp( -0.1α * d)     # probability to be absorbed while propagating a distance d [mm] *α is in [cm^-1]
    return rand() <= A ? true : false
end


# reflectivity of ray at a given plane (isotropic)
function get_reflectivity(ray::Ray, crystal::Crystal_isotropic, normal::Vector{Float64})::Float64
    θi = acos(dot(-normal, ray.k) / (norm(normal) * norm(ray.k)))
    if θi > asin(1 / crystal.n)
        return 1.0
    end
    θo = asin(crystal.n * sin(θi))      # Snell's law
    m = cross(ray.k, normal)            # vector normal to incident plane
    θpol = acos(abs(dot(ray.E, m )))    # angle between E and m
    ratio_Es = cos(θpol) / (cos(θpol) + sin(θpol))
    ratio_Ep = 1 - ratio_Es

    ts = 2 * sin(θo) * cos(θi) / sin(θi + θo)
    tp = 2 * sin(θo) * cos(θi) / (sin(θi + θo) * cos(θi - θo))
    T = tan(θi) / tan(θo) * (ratio_Es * ts^2 + ratio_Ep * tp^2)
    return 1 - T
end


# reflectivity of ray at a given plane (uniaxial)
function get_reflectivity(ray::Ray, crystal::Crystal_uniaxial, normal::Vector{Float64})::Float64
    θ = acos( abs( dot(ray.E, crystal.caxis) ) )    # angle of polarization with respect to c-axis
    ratio_e = cos(θ) / (cos(θ) + sin(θ))
    ratio_o = 1 - ratio_e

    θi = acos(dot(-normal, ray.k) / (norm(normal) * norm(ray.k)))
    # ordinary ray component (perpendicular to c-axis)
    if θi > asin(1 / crystal.no)
        Ro = 1.0
    else
        θo = asin(crystal.no * sin(θi))      # Snell's law
        m = cross(ray.k, normal)            # vector normal to incident plane
        θpol = acos(abs(dot(ray.E, m)))    # angle between E and m
        ratio_Es = cos(θpol) / (cos(θpol) + sin(θpol))
        ratio_Ep = 1 - ratio_Es
        ts = 2 * sin(θo) * cos(θi) / sin(θi + θo)
        tp = 2 * sin(θo) * cos(θi) / (sin(θi + θo) * cos(θi - θo))
        Ro = 1 - tan(θi) / tan(θo) * (ratio_Es * ts^2 + ratio_Ep * tp^2)
    end
    # extraordinary ray component (parallel to c-axis)
    if θi > asin(1 / crystal.ne)
        Re = 1.0
    else
        θo = asin(crystal.ne * sin(θi))      # Snell's law
        m = cross(ray.k, normal)            # vector normal to incident plane
        θpol = acos(abs(dot(ray.E, m)))    # angle between E and m
        ratio_Es = cos(θpol) / (cos(θpol) + sin(θpol))
        ratio_Ep = 1 - ratio_Es
        ts = 2 * sin(θo) * cos(θi) / sin(θi + θo)
        tp = 2 * sin(θo) * cos(θi) / (sin(θi + θo) * cos(θi - θo))
        Re = 1 - tan(θi) / tan(θo) * (ratio_Es * ts^2 + ratio_Ep * tp^2)
    end
    return Ro * ratio_o + Re * ratio_e
end


# boolean function to judge if ray is reflected
function judge_reflection(ray::Ray, crystal, normal::Vector{Float64})::Bool
    R = get_reflectivity(ray, crystal, normal)
    return rand() <= R ? true : false
end


# ray tracing
function ray_tracing!(ray::Ray, crystal, Δd::Float64, max_count)::Bool
    while ray.Nr < max_count
        flag::Bool = false
        t, normal = find_intersection(ray, crystal)
        N = div(t, Δd)
        if N != 0
            for i in 1:N
                if judge_absorbed(ray, crystal, Δd)
                    # if absorbed
                    if judge_reemit(ray, crystal)
                        redirect!(ray)
                        reset_wavelength!(ray, crystal)
                        flag = true
                        break
                    else
                        return false    # disappear
                    end
                else
                    # if not absorbed
                    forward!(ray, Δd)
                end
            end
        end
        if !flag
            if judge_absorbed(ray, crystal, t - Δd * N)
                # if absorbed
                if judge_reemit(ray, crystal)
                    redirect!(ray)
                    reset_wavelength!(ray, crystal)
                else
                    return false    # disappear
                end
            else
                # if not absorbed
                forward!(ray, t - N * Δd - SMALL)   # forward ray to a plane
                if judge_reflection(ray, crystal, normal)
                    reflect!(ray, normal)
                else
                    transmit!(ray, crystal, normal)
                    break
                end
            end
        end
    end
    return true     # ray escapes
end


# define a cosine corrector
function create_cosine_corrector(p0::Vector{Float64}, normal::Vector{Float64}, d::Float64)::CosineCorrector
    return CosineCorrector(p0, normal, d)
end


# judge if an escaping ray hits the cosine corrector and collected
function judge_detected(ray::Ray, cosine_corrector::CosineCorrector)::Bool
    t = dot((cosine_corrector.r0 - ray.r), cosine_corrector.normal) / dot(ray.k, cosine_corrector.normal)
    if t > 0
        hit_pos = ray.r + t * ray.k
        if norm(hit_pos - cosine_corrector.r0) <= (cosine_corrector.d * 0.5)
            val = abs( dot(ray.k, cosine_corrector.normal) )
            if rand() <= val
                return true
            end
        end
    end
    return false
end


# judge if the laser cooling ion reemits (internal quantum efficiency)
function judge_reemit(ray::Ray, crystal)::Bool
    αr = get_absorption_coefficient(ray, crystal) # [cm^-1]
    ηabs = αr / (αr + crystal.αb)
    if rand() <= ηabs
        if rand() <= crystal.QE
            return true
        end
    end
    return false
end;
