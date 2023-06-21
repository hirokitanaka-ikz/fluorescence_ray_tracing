struct Crystal_isotropic
    W::Float64      # width [mm]
    H::Float64      # height [mm]
    L::Float64      # length [mm]
    n::Float64      # refractive index
    θ::Float64      # corner angle [rad]
    T::Int64        # temperature [K]
    conc::Float64   # doping concentration [cm^-3]
    λ_vector::Vector{Float64}   # vector of wavelength [nm]
    If::Vector{Float64}       # vector of fluorescence intensity [a.u.]
    α::Vector{Float64}        # vector of absorption coefficient [cm^-1]
    p_planes::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}   # points (arbitrary) on each plane (order is 1:x=0, 2:x=W, 3:y=0, 4:y=H, 5:z=0, 6:z=L)
    plane_normals::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}  # normal vector of each palen
    QE::Float64     # internal quantum efficiency
    αb::Float64     # background absorption coefficient [cm^-1]
    αs::Float64     # scattering coefficient [cm^-1]
end


struct Crystal_uniaxial
    W::Float64      # width [mm]
    H::Float64      # height [mm]
    L::Float64      # length [mm]
    caxis::Vector{Float64}  # unit vector indicating the c-axis
    no::Float64     # refractive index for ordinary ray
    ne::Float64     # refractive index for extraordinary ray
    θ::Float64      # corner angle [rad]
    T::Int64        # temperature [K]
    conc::Float64   # doping concentration [cm^-3]
    λ_vector::Vector{Float64}   # vector of wavelength [nm]
    If_π::Vector{Float64}       # vector of fluorescence intensity for π-polarization [a.u.]
    If_σ::Vector{Float64}       # vector of fluorescence intensity for σ-polarization [a.u.]
    α_π::Vector{Float64}        # vector of absorption coefficient for π-polarization [cm^-1]
    α_σ::Vector{Float64}        # vector of absorption coefficient for σ-polarization [cm^-1]
    p_planes::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}   # points (arbitrary) on each plane (order is 1:x=0, 2:x=W, 3:y=0, 4:y=H, 5:z=0, 6:z=L)
    plane_normals::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}  # normal vector of each palen
    QE::Float64     # internal quantum efficiency
    αb::Float64     # background absorption coefficient [cm^-1]
    αs::Float64     # scattering coefficient [cm^-1]
end