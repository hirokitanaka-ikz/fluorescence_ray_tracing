struct Beam
    λ::Float64                  # wavelength of laser beam [nm]
    k::Vector{Float64}          # direction of laser beam
    d::Float64                  # beam diameter [mm]
    p0::Vector{Float64}         # point where beam passes
    p1::Vector{Float64}         # first intersection point
    p2::Vector{Float64}         # second intersection point
    path_length::Float64        # beam path length [mm]
    E::Vector{Float64}          # polarization of laser
    distribution::String        # beam distribution, "gauss" or "tophat"
end