using Random

physical["N"] = (32,32)
physical["β"] = 12.8

meta["meta_enabled"] = false
meta["Qlims"] = (-8,8)
meta["δq"] = 1e-3
meta["Qthr"] = (-7.999,7.999)
meta["w"] = 1e-5
meta["k"] = 1000
meta["is_static"] = [false]

sim["Ntherm"] = 1_000
sim["Nsweeps"] = 1_000
sim["initial"] = "cold"
sim["parallel_tempering"] = nothing
sim["swap_every"] = nothing

mc["ϵ_metro"] = 0.14
meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Wilson_loop_x16","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Polyakov_loop","measure_every" => 10)]
                       
system["veryverbose"] = false 
system["randomseeds"] = [Random.Xoshiro()]
system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])"
system["logfile"] = "Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"]).txt"
system["measure_dir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])"
system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])"
system["biasfile"] = "Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_nostatic"
system["usebiases"] = [
    "./metapotentials/N$(physical["N"])_beta$(physical["β"])_Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w0.0001_k1000_TEMPER1_HIGHSTAT.txt"]
