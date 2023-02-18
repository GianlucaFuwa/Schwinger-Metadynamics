module Measurements

    import ..Gaugefields: Gaugefield,plaquette,plaquette_sum
    import ..Observables: MetaCharge,TopCharge,wilson_loop_all,poly_loop_avg
    import ..Metadynamics: Bias_potential,penalty_potential

    defaultmeasures = Array{Dict,1}(undef,2)
    for i=1:length(defaultmeasures)
        defaultmeasures[i] = Dict()
    end
	
	defaultmeasures[1]["methodname"] = "Continuous_charge"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[2]["methodname"] = "Topological_charge"
    defaultmeasures[2]["measure_every"] = 1

    struct Measurement_set
        nummeas::Int64
        meas_calls::Array{Dict,1}
        meas_files::Array{IOStream,1}

        function Measurement_set(measure_dir;meas_calls=defaultmeasures,instance="")
            nummeas = length(meas_calls)
            meas_files = Array{IOStream,1}(undef,nummeas)
            for i=1:nummeas
                method = meas_calls[i]
                
                meas_overwrite = "w"
                if method["methodname"] == "Plaquette"
                    meas_files[i] = open(measure_dir*"/Plaquette"*instance*".txt",meas_overwrite)
                elseif method["methodname"] == "Action"
                    meas_files[i] = open(measure_dir*"/Action"*instance*".txt",meas_overwrite)
                elseif method["methodname"] == "Continuous_charge"
                    meas_files[i] = open(measure_dir*"/Continuous_charge"*instance*".txt",meas_overwrite)
                elseif method["methodname"] == "Topological_charge"
                    meas_files[i] = open(measure_dir*"/Topological_charge"*instance*".txt",meas_overwrite)
                elseif method["methodname"] == "Polyakov_loop"
                    meas_files[i] = open(measure_dir*"/Polyakov_loop"*instance*".txt",meas_overwrite)
                elseif occursin( "Wilson_loop", method["methodname"])
                    meas_files[i] = open(measure_dir*"/"*method["methodname"]*instance*".txt",meas_overwrite)
                else 
                    error("$(method["methodname"]) is not supported")
                end
            end
            return new(nummeas,meas_calls,meas_files)
        end
    end

    function measurements(itr,field::Gaugefield,measset::Measurement_set)
        for i = 1:measset.nummeas
            method = measset.meas_calls[i]
            measfile = measset.meas_files[i]
            if itr % method["measure_every"] == 0
                if method["methodname"] == "Plaquette"
                    plaq_re, plaq_im = plaquette_sum(field)
                    plaq_re /= field.NV
                    plaq_im /= field.NV
                    println(measfile,"$itr $plaq_re $plaq_im # plaq_re plaq_im")
                elseif method["methodname"] == "Action"
                    s = field.Sg/field.NV
                    println(measfile,"$itr $s # action")
                elseif method["methodname"] == "Continuous_charge" 
                    q = field.CV
                    println(measfile,"$itr $q # continuous charge")
                elseif method["methodname"] == "Topological_charge"
                    qt = TopCharge(field)
                    println(measfile,"$itr $qt $(qt^2) # Qtop Qtop^2")
                elseif method["methodname"] == "Polyakov_loop"
                    poly_re, poly_im = poly_loop_avg(field)
                    println(measfile,"$itr $poly_re $poly_im # poly_re poly_im")
                elseif occursin( "Wilson_loop", method["methodname"])
                    num = filter.(isdigit,method["methodname"])
                    LT = parse(Int64,num)
                    wils_re,wils_im = wilson_loop_all(field,LT)
                    wils_re = round.(wils_re,sigdigits=4)
                    wils_im = round.(wils_im,sigdigits=4)
                    println(measfile,"$itr $wils_re $wils_im # wilson_re wilson_im")
                else 
                    error("$(method["methodname"]) is not supported")
                end
                flush(measfile)
            end
        end
        return nothing
    end

    function build_measurements(itr,field::Gaugefield,buildmeasures::Measurement_set,instance::Int64)
        for i = 1:buildmeasures.nummeas
            method = buildmeasures.meas_calls[i]
            measfile = buildmeasures.meas_files[i]
            if itr % method["measure_every"] == 0
                if method["methodname"] == "Continuous_charge"
                    q = field.Q
                    println(measfile,"$itr $q # continuous charge")
                elseif method["methodname"] == "Topological_charge"
                    qt = calc_topcharge(field)
                    println(measfile,"$itr $qt $(qt^2) # Qtop TopSusc")
                end
                flush(measfile)
            end
        end
        return nothing
    end

    function calc_weights(q_vals::Array{Float64,1},b::Bias_potential)
        @inline index(q,qmin,dq) = round(Int,(q-qmin)/dq+0.5,RoundNearestTiesAway)
        weights = zeros(length(q_vals))
        for i=1:length(q_vals)
            idx = index(q_vals[i],b.Qmin,b.δq)
            weights[i] = exp(b[idx]+penalty_potential(q_vals[i],b.Qmin_thr,b.Qmax_thr,b.k))
        end
        return weights
    end
    
end