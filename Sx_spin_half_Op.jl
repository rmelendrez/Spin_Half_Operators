# Ronald Melendrez
# Sᵢˣ on site operator for s = ½
# 5/20/2022


# to construct module in julia 

module Sx_onsite_Mod
	
	struct Sx_onsite
		sites::Int64
		site::Int64
		spin::Float64
#		state::Array{ComplexF64}
		
		function Sx_onsite()
			new(1,1,1)
		end
	
		function Sx_onsite(sites::Int64,site::Int64,spin::Float64)
			allsites = UnitRange(1,sites)
			fullbasislist = UnitRange(0,(2^sites)-1)
			dim = Int64(2*spin+1)^sites
			

			# dictionary for spin ½ many body hilbert space
			fullbasislookup = Dict()
			
			for (i,b) in enumerate(fullbasislist)
        			fullbasislookup[b]=i
			end
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

			
			# <Sᵢˣ> operator
			Sx = zeros(ComplexF64,dim,dim)
			
			for element in fullbasislist
				
				L = fullbasislookup[element]
				
				transcient = spinflip(element,site-1)
				R = fullbasislookup[transcient]
				
				Sx[L,R] += 0.5
			
			end
			
			return Sx

		end

	end

end

			
		

