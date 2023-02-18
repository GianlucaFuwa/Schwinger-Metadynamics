module Local
	using Random
	using Printf
	#using Base.Threads:@threads,nthreads,threadid
	
	import ..System_parameters: Params
	import ..Gaugefields: Gaugefield,recalc_Sg!,recalc_CV!,dqar,daction
	import ..Metadynamics: Bias_potential,update_bias!,penalty_potential,index
	import ..Verbose_print: Verbose_,println_verbose
	
	function metropolis!(field::Gaugefield,ix::Int64,it::Int64,μ::Int64,rng::Xoshiro,ϵ::Float64)
		dU = randn(rng)*ϵ
		ΔSg = daction(field,ix,it,μ,dU)
		ΔCV = dqar(field,ix,it,μ,dU)
		accept = rand(rng) ≤ exp(-ΔSg)
		if accept  
			@inbounds field.U[ix,it,μ] += dU
			field.Sg += ΔSg
			field.CV += ΔCV
		end
		return accept
	end
	
	function metropolis_meta!(field::Gaugefield,bias::Bias_potential,ix::Int64,it::Int64,μ::Int64,rng::Xoshiro,ϵ::Float64,static::Bool)
		dU = randn(rng)*ϵ
		ΔCV = dqar(field,ix,it,μ,dU)
		
		old_ind = index(bias,field.CV)
		prop_ind = index(bias,field.CV+ΔCV)
		
		ΔSg = daction(field,ix,it,μ,dU) 
		ΔV = bias[prop_ind]-bias[old_ind]
		ΔV_pen = penalty_potential(bias,field.CV+ΔCV)-penalty_potential(bias,field.CV)
		accept = rand(rng) ≤ exp(-ΔSg-ΔV-ΔV_pen) 
		if accept 
			@inbounds field.U[ix,it,μ] += dU
			field.CV += ΔCV 
			field.Sg += ΔSg
			if ~static 
				update_bias!(bias,field.CV)
			end
		end
		return accept
	end

	function sweep!(field::Gaugefield,rng::Xoshiro,ϵ::Float64)
		NX,NT = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis!(field,ix,it,2,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis!(field,ix,it,1,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return numaccepts
	end

	function sweep_meta!(field::Gaugefield,bias::Bias_potential,rng::Xoshiro,ϵ::Float64,static::Bool)
		NX,NT = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis_meta!(field,bias,ix,it,2,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis_meta!(field,bias,ix,it,1,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return sum(numaccepts)
	end

end
