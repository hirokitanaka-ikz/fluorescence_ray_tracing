struct CosineCorrector
    r0::Vector{Float64}         # center point on the corrector plane
    normal::Vector{Float64}     # vector normal to the corrector plane
    d::Float64                  # diameter of the detector area
end