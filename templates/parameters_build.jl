using Random

physical["N"] = (32,32)
physical["β"] = 12.8

meta["meta_enabled"] = true
meta["abs_for_CV"] = false
meta["Qlims"] = (-7,7)
meta["δq"] = 1e-3
meta["Qthr"] = (-6.999,6.999)
meta["w"] = 1e-4
meta["k"] = 1000

sim["Ntherm"] = 50_000
sim["Nsweeps"] = 2_000_000
sim["initial"] = "cold"
sim["parallel_tempering"] = nothing
sim["swap_every"] = nothing

mc["ϵ_metro"] = 0.14
meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action","measure_every" => 10)]
                          #Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Wilson_loop_x16","measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Polyakov_loop","measure_every" => 10)]
                       
system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]
system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])"
system["logfile"] = "Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"]).txt"
system["measure_dir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_build"
system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])"
system["biasfile"] = "Qlims$(meta["Qlims"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_build"