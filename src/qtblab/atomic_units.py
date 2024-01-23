from numpy import pi

class au:
	eV= 27.211      # Hartree to eV
	kcalperMol= 627.509       # Hartree to kcal/mol
	bohr= 0.52917721    # Bohr to Angstrom
	Mprot= 1836.15        # proton mass
	hbar= 1.          # Planck's constant
	fs= 2.4188843e-2  # AU time to femtoseconds
	ps=fs/1000
	kelvin = 3.15774e5    # Hartree to Kelvin
	THz = 1000./fs		# AU frequency to THz
	nNewton = 82.387		# AU Force to nNewton 
	cm1=219471.52  # Hartree to cm-1 
	gmol_Afs=bohr/(Mprot*fs) #AU momentum to (g/mol).A/fs
	kbar=294210.2648438959 # Hartree/bohr**3 to kbar
	
	mapping={
		"1": 1.,
		"au": 1.,
		"ev": eV,
		"kcalpermol": kcalperMol,
		"angstrom": bohr,
		"bohr": 1./bohr,
		"amu": 1./Mprot,
		"femtoseconds": fs,
		"fs": fs,
		"picoseconds": fs/1000.,
		"ps":fs/1000.,
		"kelvin": kelvin,
		"k": kelvin,
		"thz": THz,
		"tradhz": THz/(2*pi),
		"cm-1":cm1,
		"kbar":kbar,
		"gmolafs":gmol_Afs,
	}

	@staticmethod
	def get_multiplier(unit_string):
		unit_string=unit_string.lower().strip()

		multiplier=1.

		unit_start=0
		unit_stop=0
		in_power=False
		power_start=0
		power_stop=0
		tmp_power=1.
		for i in range(len(unit_string)):

			current_char=unit_string[i]
			#print(current_char)
			if current_char=="^":

				if in_power:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError
					
				in_power=True
				unit_stop=i-1
				if unit_stop-unit_start < 0:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

			elif current_char=="{":

				if not in_power:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError
				
				if i+1 >= len(unit_string):
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

				power_start=i+1

			elif current_char=="}":

				if not in_power:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError
				
				in_power=False
				power_stop=i-1

				if power_stop-power_start < 0:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError
				else:
					tmp_power_read=int(unit_string[power_start:power_stop+1])
					tmp_power=tmp_power*tmp_power_read

				unit_substring=unit_string[unit_start:unit_stop+1]
				tmp_unit=au.mapping[unit_substring]

				multiplier=multiplier*(tmp_unit**tmp_power)

				power_start=0
				power_stop=0
				unit_start=i+1
				unit_stop=0
			
			elif current_char=="*":
				if in_power:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

				if unit_start==i:
					unit_start=i+1
					tmp_power=1.
					continue

				unit_stop=i-1
				if unit_stop-unit_start < 0:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

				unit_substring=unit_string[unit_start:unit_stop+1]
				tmp_unit=au.mapping[unit_substring]
				multiplier=multiplier*(tmp_unit**tmp_power)

				unit_start=i+1
				unit_stop=0
				tmp_power=1.
			
			elif current_char=='/':
				if in_power:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

				if unit_start==i:
					unit_start=i+1
					tmp_power=-1.
					continue

				unit_stop=i-1
				if unit_stop-unit_start < 0:
					print("Error: Syntax error in unit '"+unit_string+"' !")
					raise ValueError

				unit_substring=unit_string[unit_start:unit_stop+1]
				tmp_unit=au.mapping[unit_substring]
				multiplier=multiplier*(tmp_unit**tmp_power)

				unit_start=i+1
				unit_stop=0
				tmp_power=-1.

			else:
				if i+1 >= len(unit_string):
					if in_power:
						print("Error: Syntax error in unit '"+unit_string+"' !")
						raise ValueError

					unit_stop=i
					if unit_stop-unit_start < 0 :
						print("Error: Syntax error in unit '"+unit_string+"' !")
						raise ValueError

					unit_substring=unit_string[unit_start:unit_stop+1]
					tmp_unit=au.mapping[unit_substring]
					multiplier=multiplier*(tmp_unit**tmp_power)
			
		return multiplier
