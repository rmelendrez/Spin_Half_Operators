# Ronald Melendrez
# Sʸᵢ on site operator for s = ½
# 5/20/2022


# to construct module in julia 

module Sy_onsite_Mod
	
	struct Sy_onsite
		
		# struct parameters
		sites::Int64
		site::Int64
		spin::Float64
		
		# initiate function
		function Sy_onsite()
			new(1,1,1)
		end
	
		function Sy_onsite(sites::Int64,site::Int64,spin::Float64)
			# construct sites
			allsites = UnitRange(1,sites)
			# construct basis list
			fullbasislist = UnitRange(0,(2^sites)-1)
			# dimension of hilbert space for spin ½
			dim = Int64(2*spin+1)^sites
			

			# initiate dictionary for spin ½ many body hilbert space
			fullbasislookup = Dict()
			
			# construct library for basis look up
			for (i,b) in enumerate(fullbasislist)
				fullbasislookup[b] = i
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

			
			# < Sᶻᵢ > operator
			Sy = zeros(ComplexF64,dim,dim)
			
			# Goes through the entire hilbert space
			for element in fullbasislist
				
				# Gets the column index for the matrix
				#  <R |  |R'> Sʸᵢ <L'|  |L>	
				L = fullbasislookup[element]
				
				# Sʸᵢ = (S⁺ᵢ - S⁻ᵢ)/(2i)
				if getspin(element,site-1) == 1
					sign = -1.0
				else
					sign = +1.0
				end
				
				transcient = spinflip(element,site-1)
				R = fullbasislookup[transcient]
				
				Sy[R,L] += 0.5*(-im)*sign
			
			end
			
			return Sy

		end

	end

end

			
		

