module Mainbuild
    using Printf 
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    using Base.Threads
    
    import ..System_parameters: Params,Params_set,parameterloading
    import ..Verbose_print: Verbose_,println_verbose,print2file
    import ..Gaugefields: Gaugefield,recalc_Sg!,recalc_CV!
    import ..Metadynamics: Bias_potential
    import ..Measurements: Measurement_set,measurements,calc_weights
    import ..Local: sweep!,sweep_meta!
    import ..Tempering: tempering_swap!

    import ..System_parameters:physical,meta,sim,mc,meas,system

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,sim,mc,meas,system)

        run_build(params_set)

        return nothing
    end

    function run_build(params_set::Params_set)
        params = parameterloading(params_set)

        field = Gaugefield(params)
        bias = Bias_potential(params)
        run_build!(field,bias,params)

        return nothing
    end

    function run_build!(field::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)

        measset = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        rng = params.randomseeds[1]

        if params.initial == "hot"
            field.U = rand(size(field.U) .- 0.5)*2*2pi
            recalc_Sg!(field)
            recalc_CV!(field)
        end 

        for therm = 1:params.Ntherm
            sweep!(field,rng,ϵ)
        end
        
        numaccepts = 0
        for trj = 1:params.Nsweeps
            tmp = sweep_meta!(field,bias,rng,ϵ,false)
            numaccepts += tmp
            # writing logs....#
            if params.veryverbose
            println_verbose(verbose," ",trj," ",100*numaccepts/trj/field.NV/2,"%"," # trj accrate")
            end
            #-----------------#
            measurements(trj,field,measset)
        end
        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.NV/2)%")

        writedlm(bias.fp,[bias.q_vals bias.values])
        flush(bias)

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
