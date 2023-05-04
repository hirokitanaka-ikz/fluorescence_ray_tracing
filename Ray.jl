mutable struct Ray
    r::Vector{Float64}      # position: r = [x, y, z] [mm]
    r0::Vector{Float64}     # initial position: r0 = [x0, y0, z0] [mm]
    k::Vector{Float64}      # Poynthing vector: k = [kx, ky, kz]
    E::Vector{Float64}      # Polarization vector: E = [Ex, Ey, Ez]
    λ::Float64              # current wavelength [nm]
    λgen::Float64           # generated wavelength [nm]
    Nr::Int64               # number of reflections
    Ray(r, k, E, λ) = new(r, r, k, E, λ, λ, 0)
end
