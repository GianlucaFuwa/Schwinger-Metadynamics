module Metadynamics
	using DelimitedFiles
	import ..System_parameters: Params
	import ..Gaugefields: Gaugefield,plaquette
	
	mutable struct Bias_potential
		abs_for_CV::Bool
		Qlims::NTuple{2,Float64}
		Qthr::NTuple{2,Float64}
		δq::Float64
		w::Float64
		k::Int64
		is_static::Bool
		values::Array{Float64,1}
		q_vals::Array{Float64,1}
		fp::IOStream
		
		function Bias_potential(p::Params;instance::Int64=1)
			abs_for_CV = p.abs_for_CV
			Qlims = p.Qlims
			Qthr = p.Qthr
			δq = p.δq
			w = p.w
			k = p.k
			is_static = p.is_static[instance]
			values = potential_from_file(p,p.usebiases[instance])
			q_vals = range(p.Qlims[1],p.Qlims[2],Int((p.Qlims[2]-p.Qlims[1])/p.δq)+1)
			fp = open(p.biasfiles[instance],"w")
			return new(abs_for_CV,Qlims,Qthr,δq,w,k,is_static,values,q_vals,fp)
		end	
	end

	function potential_from_file(p::Params,usebias::Union{Nothing,String})
		if usebias === nothing
			println("usebias is nothing")
			return zeros(round(Int,(p.Qlims[2]-p.Qlims[1])/p.δq,RoundNearestTiesAway)+1)
		else
			values = readdlm(usebias,Float64)
			@assert length(values[:,2]) == round(Int,((p.Qlims[2]-p.Qlims[1])/p.δq),RoundNearestTiesAway)+1 "Potential doesn't fit its parameters"
			return values[:,2]
		end
	end

	function Base.flush(b::Bias_potential)
        if b.fp !== nothing
            flush(b.fp)
        end
    end

	function Base.seekstart(b::Bias_potential)
		if b.fp !== nothing
            seekstart(b.fp)
        end
    end

	function Base.setindex!(b::Bias_potential,v,i::Int)
		b.values[min(length(b.values),max(i,1))] = v
		return nothing
	end

	@inline function Base.getindex(b::Bias_potential,i::Int)
		return b.values[min(length(b.values),max(i,1))]
	end

	@inline function index(b::Bias_potential,q::Float64)
		q = b.abs_for_CV ? abs(q) : q
		grid_index = (q-b.Qlims[1])/b.δq + 0.5
		return round(Int,grid_index,RoundNearestTiesAway)
	end

	function update_bias!(b::Bias_potential,q::Float64)
		q = b.abs_for_CV ? abs(q) : q
		grid_index = index(b,q)
		grid_q = b.Qlims[1] + grid_index*b.δq

		b[grid_index] += b.w*(1-(q-grid_q)/b.δq)
		b[grid_index+1] += b.w*(q-grid_q)/b.δq
		return nothing
	end

	function penalty_potential(b::Bias_potential,q::Float64)
		q = b.abs_for_CV ? abs(q) : q
		if q < b.Qthr[1] || q > b.Qthr[2]
			p_pot = b.k*min((q-b.Qthr[1])^2,(q-b.Qthr[2])^2)
		else 
			p_pot = 0
		end
		return p_pot
	end

end


