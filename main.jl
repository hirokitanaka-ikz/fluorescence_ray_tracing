using Printf
using Distributed
using CSV, DataFrames
using Dates
using YAML
using OrderedCollections
using Statistics

include("simulation.jl")


function write_log(logfile, message, with_timestamp::Bool=false)
    if with_timestamp
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        full_message = "[$timestamp] $message\n"
    else
        full_message = "$message\n"
    end
    # output to terminal
    print(full_message)

    # output to file
    open(logfile, "a") do f
        write(f, full_message)
    end
end


function print_parameters(logfile, params)
    indent = "  "
    # print to Terminal
    write_log(logfile, "-----Crystal-----")
    for (key, value) in params["crystal"]
        write_log(logfile, indent * "$key: $value")
    end
    write_log(logfile, "\n-------Beam-------")
    for (key, value) in params["beam"]
        write_log(logfile, indent * "$key: $value")
    end
    write_log(logfile, "\n-----Simulation-----")
    for (key, value) in params["simulation"]
        write_log(logfile, indent * "$key: $value")
    end
    write_log(logfile, "*****(Parameters)*****")
end


function main()

    # create directory
    if !isdir("result")
        mkdir("result")
    end
    dt = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkdir("./result/$dt")
    logfile = "./result/$dt/log.txt"

    # print first message
    first_message = "-------------------------------------------------------------\n\
    | Fluorescence Ray Tracing Simulation by Monte Carlo Method |\n\
    -------------------------------------------------------------\n"
    print(first_message)
    open(logfile, "a") do f
        write(f, first_message)
    end

    # load parameters from yaml
    params = YAML.load_file("params.yaml", dicttype=OrderedDict{String, Any})
    write_log(logfile, "*****Parameters*****")
    print_parameters(logfile, params)

    # run simulation
    nproc = length(procs())
    write_log(logfile, "Simulation Start ($nproc Processes)", true)
    elapsed = @elapsed df = simulation(params)
    write_log(logfile, "Simulation Successfully Finish (Elapsed Time: $elapsed sec.)", true)
    N = convert(Int64, params["simulation"]["N"])
    Nesc = size(df, 1)

    write_log(logfile, "\n----------Result----------")
    write_log(logfile, "Escaping Fluorescence Rays: $Nesc / $N\nFluorescence Escape Efficiency = $(Nesc / N * 100)%")
    mean_wl_gen = mean(df.λgen)
    std_wl_gen = std(df.λgen)
    mean_wl_esc = mean(df.λesc)
    std_wl_esc = std(df.λesc)
    write_log(logfile, "Mean Fluorescence Wavelength (Generate): $mean_wl_gen +- $std_wl_gen nm")
    write_log(logfile, "Mean Fluorescence Wavelength (Escape): $mean_wl_esc +- $std_wl_esc nm")


    # save csv (result)
    rename!(df, :λgen => :wavelength_gen)
    rename!(df, :λesc => :wavelength_esc)
    csv_filepath = "result/$dt/out.csv"
    CSV.write(csv_filepath, df)

    # save yaml (parameters)
    cp("params.yaml", "result/$dt/params.yaml")
    println("\nPROGRAM FINISH")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end