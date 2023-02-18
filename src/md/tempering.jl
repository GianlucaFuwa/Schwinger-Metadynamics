module Tempering
    using Random
    
    import ..Gaugefields: Gaugefield
    import ..Metadynamics: Bias_potential,index,update_bias!,penalty_potential

    function tempering_swap!(field_main::Gaugefield,field_meta::Gaugefield,bias::Bias_potential,rng::Xoshiro,static::Bool)
		main_ind = index(bias,field_main.CV)
		meta_ind = index(bias,field_meta.CV)

		ΔCV = field_meta.CV - field_main.CV
		ΔSg = field_meta.Sg - field_main.Sg
		ΔV = bias[meta_ind] - bias[main_ind]
		ΔV_pen = penalty_potential(bias,field_meta.CV) - penalty_potential(bias,field_main.CV)
		accept_swap = rand(rng) ≤ exp(ΔV+ΔV_pen)
		if accept_swap
			swap!(field_main,field_meta)
			field_meta.CV -= ΔCV
			field_main.CV += ΔCV
			field_meta.Sg -= ΔSg
			field_main.Sg += ΔSg
			if ~static
				update_bias!(bias,field_meta.CV)
			end
		end
		return accept_swap
	end 

end