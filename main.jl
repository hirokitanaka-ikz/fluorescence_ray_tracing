using Printf
using Distributed
using CSV, DataFrames
using Dates
using YAML
# using FileIO

include("simulation.jl")
# @everywhere include("parameters.jl")

function main()
    # load parameters from yaml
    params = YAML.load_file("params.yaml")
    elapsed = @elapsed df = simulation(params)
    dt = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    if !isdir("result")
        mkdir("result")
    end
    mkdir("./result/$dt")

    FILENAME = "$(params["crystal"]["name"])_$(params["crystal"]["doping_level"])%_" *
    "$(params["crystal"]["W"])x$(params["crystal"]["H"])x$(params["crystal"]["L"])_" *
    "$(params["crystal"]["shape"])_$(params["beam"]["distribution"])"

    # save csv (result)
    csv_filepath = "result/$dt/$FILENAME.csv"
    CSV.write(csv_filepath, df)

    # save yaml (parameters)
    cp("params.yaml", "result/$dt/params.yaml")
    @printf("Elapsed time: %.3f sec.\n", elapsed)
end