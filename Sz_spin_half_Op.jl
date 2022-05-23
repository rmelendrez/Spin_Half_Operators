# Ronald Melendrez
# Sᵢᶻ on site operator for s = ½
# 5/20/2022


# to construct module in julia 

module Sz_onsite_Mod
	
	struct Sz_onsite
		sites::Int64
		site::Int64
		spin::Float64
		state::Array{ComplexF64}
		
		function Sz_onsite()
			new(1,1,1)
		end
	
		function Sz_onsite(sites::Int64,site::Int64,spin::Float64)
			allsites = UnitRange(1,sites)
			fullbasislist = UnitRange(0,(2^sites)-1)
			dim = Int64(2*spin+1)^sites
			

			# dictionary for spin ½ many body hilbert space
			fullbasislookup = Dict()
						
			# Sᵢᶻ get spin projection along z in binary notation
			# 1 = |↑>  and  0 = |↓>
			#  
			function getspin(b,i)
				b >> i & 1
			end
		
			# Sᵢ⁺ and Sᵢ⁻ spin flip operator
			function spinflip(b,i)
				if ( b >> i ) & 1 == 1
					return b - ( 1 << i )
				else
					return b + ( 1 << i )
				end
			end

			for (i,b) in enumerate(fullbasislist)
        			fullbasislookup[b]=i
			end
			# <Sᵢˣ> operator
			Sz = zeros(ComplexF64,dim,dim)
			
			for element in fullbasislist
				
			L = fullbasislookup[element]
                        
				if getspin(element,site-1) == 1
                                	sign = 1.0
                        	else
                                	sign = -1.0
                        	end

				Sz[L,L] = sign*0.5
			
			end
			
			return Sz

		end

	end

end

			
		

