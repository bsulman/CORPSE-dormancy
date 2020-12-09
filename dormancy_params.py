params={
	'vmaxref':[600,25,600], #Relative maximum enzymatic decomp rates
	'Ea':[37e3,54e3,50e3],	# Activation energy
	'kC':[1.0e-6,1.0e-6,1.0e-6], # Michaelis-Menten parameter
#	'kC':[5.3e-3,5.3e-3,5.3e-3],	# Michaelis-Menton parameter
	'gas_diffusion_exp':2.5,  # Determines suppression of decomp at high soil moisture
	'minMicrobeC':0.0,	   #Minimum microbial biomass (fraction of total C)
	#'Tmic_30C':0.005,	 # 0.03  # Microbial lifetime  (years)
	'Q10_Tmic':2.0,
	'Tmic_20C':0.02, # 0.04  # Microbial lifetime  (years) #
	'Tmic_dormant':0.5, # 3.0 Dormant microbe lifetime (years)
	'et':0.3,		  # Fraction of turnover not converted to CO2
	'et_dormant':0.9,   # Fraction of dormant microbe turnover not converted to CO2
	'eup':[0.4,0.02,0.4], # Carbon uptake efficiency
	'tProtected':75.0,	# Protected C turnover time (years)
	'protection_rate':[1.0,0.0,1.0], # Protected carbon formation rate (year-1)
	'max_activation_rate':0.21*365*24 ,	# Dormant microbe max activation rate (fraction per year)
	'k_activation_20C':0.025   , #0.0003			 # Half-saturation of growth/turnover ratio in activation
	'k_activation_30C':0.008,   # 0.0003  activation of dormant biomass is faster in heated soils than in unheated soils
	'max_dormancy_rate':0.5*365*24  ,  # Active microbe max dormancy rate (fraction per year)
	'max_active_frac_30C':0.40 , #0.2 #0.15
	'max_active_frac_20C':0.40 , #0.4 #0.30, #  Max % active biomass in dry, unheated (T=20C) is higher than in dry, heated (T=30C) soils because of lesser warming-induced stress
	'min_active_frac_30C':0.0,
	'min_active_frac_20C':0.0,
	'fast_dormancy':True
}

params_nodormant=params.copy()
params_nodormant['max_dormancy_rate']=0.0;params_nodormant['max_activation_rate']=0.0;params_nodormant['k_activation_30C']=0.0;params_nodormant['k_activation_20C']=0.0
params_nodormant['max_active_frac_30C']=1.0;params_nodormant['max_active_frac_20C']=1.0
params_nodormant['et']=0.7
params_nodormant['kC']=[2.5e-4,2.5e-4,2.5e-4];params_nodormant['Tmic_20C']=0.2;#params_nodormant['Tmic_30C']=0.1

# Spun up cohort
import CORPSE
#original
#c_spunup=CORPSE.soil_carbon_cohort(litterC=[0.0094,1.7,0.0035], protectedC=[0.13,0,0.45], activeMicrobeC=0.002, dormantMicrobeC=0.0056)

# Initial biomass
# T=30 C
c_spunup=CORPSE.soil_carbon_cohort(litterC=[0.0094,1.7,0.0035], protectedC=[0.13,0,0.45], activeMicrobeC=0.00005, dormantMicrobeC=0.00455)

# T=20 C
#c_spunup=CORPSE.soil_carbon_cohort(litterC=[0.0094,1.7,0.0035], protectedC=[0.13,0,0.45], activeMicrobeC=0.00030, dormantMicrobeC=0.0072)
