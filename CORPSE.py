# This should be a complete list of parameters with descriptions
# The set_params function will check that all are present and there are no extras
expected_params={	'vmaxref': 'Relative maximum enzymatic decomp rates (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
        	'minMicrobeC':	   'Minimum microbial biomass (fraction of total C)',
        	# 'Tmic_30C':	   'Microbial lifetime at 30C (years)',
            'Q10_Tmic':  'Q10 of microbial lifetime',
        	'Tmic_20C': 'Microbial lifetime at 20C (years)',
        	'Tmic_dormant': 'Dormant microbe lifetime (years)',
        	'et':		  'Fraction of turnover not converted to CO2',
        	'et_dormant':   'Fraction of dormant microbe turnover not converted to CO2',
        	'eup': 'Carbon uptake efficiency (length 3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
        	'max_activation_rate': 'Dormant microbe max activation rate (fraction per year)',
        	'k_activation_30C': 'Half-saturation of growth/turnover ratio in activation at 30C',
        	'k_activation_20C': 'Half-saturation of growth/turnover ratio in activation at 20C',
        	'max_dormancy_rate': 'Active microbe max dormancy rate (fraction per year)',
        	'max_active_frac_30C': 'Maximum active fraction at 30C',
        	'max_active_frac_20C': 'Maximum active fraction at 20C',
        	'min_active_frac_30C': 'Minimum active fraction at 30C',
        	'min_active_frac_20C': 'Minimum active fraction at 20C',
        	'fast_dormancy': 'Whether to use instant dormancy fraction calculation (True or False)',
            }

class soil_carbon_cohort:
    def __init__(self,
            litterC=[0.0,0.0,0.0],
            protectedC=[0.0,0.0,0.0],
            activeMicrobeC=0.0,dormantMicrobeC=0.0,
            CO2=0.0,params=None):
        '''Initialize a soil carbon cohort.
           litterC: Unprotected C, available for decomposition. (kgC/m2)
           protectedC: Physically or chemically protected, inaccessible to microbes. (kgC/m2)
           activeMicrobeC: Microbes that decompose the carbon. (kgC/m2)
           dormantMicrobeC: Dormant microbes with capacity to wake up (kgC/m2)
           CO2: Cumulative CO2 production (kgC/m2)
           originalC: Cumulative carbon added to cohort, used for conservation checks. (kgC/m2)
           params: Parameter set (see set_params method)
           '''

        from numpy import array

        self.litterC=array(litterC)
        self.protectedC=array(protectedC)
        self.activeMicrobeC=activeMicrobeC
        self.dormantMicrobeC=dormantMicrobeC
        self.CO2=CO2
        self.originalC=sum(litterC+protectedC)+activeMicrobeC+dormantMicrobeC+CO2


        if params is not None:
            self.set_params(params)

        return None


    def set_params(self,params):
        '''params: dictionary containing parameter values. Should contain these fields (showing reasonable default values):
                 vmaxref=[2500,600,2000]; Relative maximum enzymatic decomp rates
                 Ea=[37e3,54e3,50e3];     Activation energy
                 kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
                 gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
                 minMicrobeC=1e-3;       Minimum microbial biomass (fraction of total C)
                 Tmic=0.15;          Microbial turnover rate
                 Tmic_dormant = 5.0; Dormant microbe turnover rate
                 et=0.5;           Fraction of turnover not converted to CO2
                 eup=[0.6,0.05,0.6];  Carbon uptake efficiency
                 tProtected=75.0;     Protected C turnover time (years)
                 protection_rate=[1.0,0.0,1.0];  Protected carbon formation rate (year-1)'''

        from numpy import iterable,array
        self.params=params
        unused_params=expected_params.copy()
        for k in self.params.keys():
            if k not in expected_params:
                raise ValueError('Parameter set contains unexpected parameter %s'%k)
            unused_params.pop(k)
            if iterable(self.params[k]):
                self.params[k]=array(self.params[k])
        if len(unused_params)>0:
            for k in unused_params.keys():
                print ('Missing parameter: %s [%s]'%(k,unused_params[k]))
            raise ValueError('Missing parameters: %s'%unused_params.keys())


    def add_carbon(self,litterC=[0.0,0.0,0.0],protectedC=[0.0,0.0,0.0],
                        activeMicrobeC=0.0,dormantMicrobeC=0.0, CO2=0.0):
        'Add carbon to the cohort, keeping track of total for conservation checks.'

        from numpy import array

        if(any(array(litterC)<0.0)):
            raise ValueError('litterC must be >= 0')
        if(any(array(protectedC)<0.0)):
            raise ValueError('protectedC must be >= 0')
        if activeMicrobeC<0.0:
            raise ValueError('activeMicrobeC must be >=0.0')
        if dormantMicrobeC<0.0:
            raise ValueError('dormantMicrobeC must be >=0.0')
        if CO2<0.0:
            raise ValueError('CO2 must be >=0.0')

        self.litterC=self.litterC+array(litterC)
        self.protectedC=self.protectedC+array(protectedC)
        self.activeMicrobeC=self.activeMicrobeC+activeMicrobeC
        self.dormantMicrobeC=self.dormantMicrobeC+dormantMicrobeC
        self.CO2=self.CO2+CO2
        self.originalC=self.originalC+sum(litterC+protectedC)+activeMicrobeC+dormantMicrobeC+CO2

    def update(self,T,theta,dt):
        '''Update the cohort, with decomposition, microbial growth, etc.
           T: Temperature (K)
           theta: Soil water content (fraction of saturation)
           dt: Time step (years)

           Returns a dictionary of outputs.'''

        if T<0.0:
            raise ValueError('T is in K and must be >=0')

        if theta<0.0:
            theta=0.0
        elif theta>1.0:
            theta=1.0

        totalResp=0.0;

        et=self.params['et']
        eup=self.params['eup']

        # Active/dormant transfers

        if not self.params['fast_dormancy']:
            activation = self.activation_rate(T,theta)
            deactivation = self.dormancy_rate(T,theta)
            self.activeMicrobeC=self.activeMicrobeC+(activation-deactivation)*dt
            self.dormantMicrobeC=self.dormantMicrobeC+(deactivation-activation)*dt
        else:
            dormantfrac=self.steady_state_dormant_frac(T,theta)

            totalMic=self.activeMicrobeC+self.dormantMicrobeC
            self.activeMicrobeC=totalMic*(1.0-dormantfrac)
#            print(self.activeMicrobeC)
            self.dormantMicrobeC=totalMic*(dormantfrac)

        # Calculate maximum potential C decomposition rate
        tempResp=self.Resp(self.litterC,self.activeMicrobeC,T,theta)

        # Carbon loss cannot exceed size of pool
        tempResp[dt*tempResp > self.litterC] = self.litterC[dt*tempResp > self.litterC]/dt


        # Microbial turnover
        Tmic=self.Tmic_T(T)
        microbeTurnover=max(0.0,(self.activeMicrobeC-self.params['minMicrobeC']*sum(self.litterC))/Tmic);   # kg/m2/yr
        dormantMicrobeTurnover=self.dormantMicrobeC/self.params['Tmic_dormant']

        maintenance_resp=microbeTurnover*(1.0-et)+dormantMicrobeTurnover*(1.0-self.params['et_dormant'])

        deadmic_C_produced=dt*(microbeTurnover*et+dormantMicrobeTurnover*self.params['et_dormant'])   # actual fraction of microbial turnover
        self.litterC[2]=self.litterC[2]+deadmic_C_produced  # kg/m2

        # CO2 production and cumulative CO2 produced by cohort
        CO2prod=dt*(sum(tempResp*(1.0-eup))+maintenance_resp) # kg/m2
        self.CO2=self.CO2+CO2prod  # kg/m2

        microbeGrowth=sum(tempResp*eup);

        self.activeMicrobeC=self.activeMicrobeC + dt*(microbeGrowth-microbeTurnover);
        self.dormantMicrobeC=self.dormantMicrobeC - dt*dormantMicrobeTurnover


        # Update the amount of organic C and N in the cohort after the decomposition process

        self.litterC=self.litterC-dt*tempResp;     # kg/m2
        totalResp=totalResp+tempResp;             # kg/m2/yr


        # Update protected carbon
        protectedCturnover = self.protectedC/self.params['tProtected'] ;

        if(sum(self.litterC)>0.0):

            newProtectedC = self.params['protection_rate']*self.litterC*dt;
            #  kg/m2      =   yr-1             kg/m2        yr

        else:
            newProtectedC = [0.0,0.0,0.0];

        newProtectedC[newProtectedC>self.litterC] = self.litterC[newProtectedC>self.litterC];
        self.protectedC = self.protectedC + newProtectedC - dt*protectedCturnover;
        self.litterC = self.litterC - newProtectedC + dt*protectedCturnover;

        protected_produced=newProtectedC;  # kg/m2
        protected_turnover_rate=protectedCturnover;  # kg/m2/dt

        outputs={}

        outputs['decomp']=totalResp;
        outputs['protected_produced']=protected_produced;
        outputs['protected_turnover_rate']=protected_turnover_rate;
        outputs['CO2prod']=CO2prod
        # outputs['activation']=activation
        # outputs['dormancy']=deactivation

        return outputs

    # Activation of dormant microbes
    def activation_rate(self,T,theta):
        potential_growth_rate=sum(self.Resp(self.litterC,sum(self.litterC*10),T,theta)*self.params['eup'])
        potential_turnover_rate=(self.activeMicrobeC+self.dormantMicrobeC)/self.params['Tmic']

        if potential_turnover_rate == 0.0:
            return 0.0

        g_t_ratio=potential_growth_rate/potential_turnover_rate
        return self.params['max_activation_rate']*g_t_ratio/(g_t_ratio+self.params['k_activation'])*self.dormantMicrobeC

    def dormancy_rate(self,T,theta):
        potential_growth_rate=sum(self.Resp(self.litterC,sum(self.litterC*10),T,theta)*self.params['eup'])
        potential_turnover_rate=(self.activeMicrobeC+self.dormantMicrobeC)/self.params['Tmic']

        if potential_growth_rate==0:
            return self.params['max_dormancy_rate']*self.activeMicrobeC

        else:

            return min(1.0,potential_turnover_rate/potential_growth_rate)*self.activeMicrobeC*self.params['max_dormancy_rate']

    def g_t_ratio(self,T,theta):
        potential_growth_rate=sum(self.Resp(self.litterC,sum(self.litterC*1),T,theta)*self.params['eup'])
        Tmic=self.Tmic_T(T)
        potential_turnover_rate=(sum(self.litterC)*1)/Tmic
        g_t_ratio=potential_growth_rate/potential_turnover_rate
        return g_t_ratio


    def steady_state_dormant_frac(self,T,theta):
        potential_growth_rate=sum(self.Resp(self.litterC,sum(self.litterC*1),T,theta)*self.params['eup'])
        # potential_turnover_rate=(self.activeMicrobeC+self.dormantMicrobeC)/self.params['Tmic']
        Tmic=self.Tmic_T(T)

        potential_turnover_rate=(sum(self.litterC)*1)/Tmic

        g_t_ratio=potential_growth_rate/potential_turnover_rate
#        print("g_t_ratio",g_t_ratio)
#        print("potential_turnover_rate",potential_turnover_rate)
#        print("potential_growth_rate",potential_growth_rate)

        max_active_frac=self.max_active_frac_T(T)
        k_activation=self.k_activation_T(T)


#        result = 1.0-max_active_frac*(g_t_ratio/(g_t_ratio+k_activation))
#        if (1-result) < self.params['min_active_frac_30C']:
#            if T==293.15:
#                result = 1-self.params['min_active_frac_20C']
#            elif T==303.15:
#                result = 1-self.params['min_active_frac_30C']

#        return result

        return 1.0-max_active_frac*(g_t_ratio/(g_t_ratio+k_activation))

    def Tmic_T(self,T):
        return self.params['Tmic_20C']/self.params['Q10_Tmic']**((T-293.15)/10.0)
        # return (self.params['Tmic_20C']+(T-293.15)*(self.params['Tmic_30C']-self.params['Tmic_20C'])/10.0)
    def max_active_frac_T(self,T):
        return self.params['max_active_frac_20C']+(T-293.15)*(self.params['max_active_frac_30C']-self.params['max_active_frac_20C'])/10.0
    def k_activation_T(self,T):
        return self.params['k_activation_20C']+(max(T,290.0)-293.15)*(self.params['k_activation_30C']-self.params['k_activation_20C'])/10.0

    # Decomposition rate
    def Resp(self,Ctotal,Chet,T,theta):
        '''Chet        heterotrophic (microbial) C, living microbial biomass in the cohort
        T,theta     temperature (k), theta (fraction of 1.0)
        Ctotal      Substrate C (3-value vector)'''


        if(sum(Ctotal)==0.0 or theta==0.0 or Chet==0.0):
            return Ctotal*0.0


        Resp=self.Vmax(T)*theta**3.0*(Ctotal)*Chet/(sum(Ctotal)*self.params['kC']+Chet)*(1.0-theta)**self.params['gas_diffusion_exp'];

        return Resp



    def Vmax(self,T):
        Tref=293.15;
        Rugas=8.314472;

        from numpy import exp

        # Normalization value
        alpha=self.params['vmaxref']/exp(-self.params['Ea']/(Rugas*Tref));
        Vmax=alpha*exp(-self.params['Ea']/(Rugas*T));
        return Vmax


    def totalC(self):
        return sum(self.litterC+self.protectedC)+self.activeMicrobeC+self.dormantMicrobeC+self.CO2;

    def check_validity(self,tolerance=1e-6):
        'Check for valid values and carbon conservation'

        from numpy import isfinite

        bad=False;
        if any(self.litterC<0) or any(self.protectedC<0):
            bad=True;

        if self.activeMicrobeC<0 or self.dormantMicrobeC<0 or self.CO2<0 or self.originalC<0:
            bad=True;

        if any(~isfinite(self.litterC)) or any(~isfinite(self.protectedC)) \
                or ~isfinite(self.activeMicrobeC) or ~isfinite(self.dormantMicrobeC) \
                or ~isfinite(self.CO2) or ~isfinite(self.originalC):
            bad=True;

        if bad:
            print(self)
            raise RuntimeError('Cohort bad')

        # Check carbon conservation
        if abs(self.totalC()-self.originalC) > tolerance:
            print(self)
            raise RuntimeError('Cohort C conservation violated. Sum = %1.2d, Original = %1.2d'%(self.totalC(),self.originalC))



    def __mul__(self,value):
        # if not is_numlike(value):
        #     raise ValueError('Must multiply by a number')

        return soil_carbon_cohort(
            litterC=self.litterC*value,
            protectedC=self.protectedC*value,
            activeMicrobeC=self.activeMicrobeC*value,
            dormantMicrobeC=self.dormantMicrobeC*value,
            CO2=self.CO2*value)

    def __rmul__(self,value):
        return(self*value)

    def __div__(self,value):
        return self*(1.0/value)

    def __add__(self,cohort2):
        if not isinstance(cohort2,soil_carbon_cohort):
            raise ValueError('Cohorts can only be added to other cohorts')

        return soil_carbon_cohort(
            litterC=self.litterC+cohort2.litterC,
            protectedC=self.protectedC+cohort2.protectedC,
            activeMicrobeC=self.activeMicrobeC+cohort2.activeMicrobeC,
            dormantMicrobeC=self.dormantMicrobeC+cohort2.dormantMicrobeC,
            CO2=self.CO2+cohort2.CO2)

    def __neg__(self):
        return self*-1

    def __sub__(self,cohort2):
        return self + -cohort2

    def copy(self):
        out=self*1.0
        if hasattr(self,'params'):
            out.params=self.params
        return out


    def __str__(self):
        s='soil_carbon_cohort(litterC=[%1.2g,%1.2g,%1.2g], protectedC=[%1.2g,%1.2g,%1.2g], activeMicrobeC=%1.2g, dormantMicrobeC=%1.2g, CO2=%1.2g, originalC=%1.2g)'\
                    %(self.litterC[0],self.litterC[1],self.litterC[2],
                      self.protectedC[0],self.protectedC[1],self.protectedC[2],
                      self.activeMicrobeC,self.dormantMicrobeC,self.CO2,self.originalC)
        return s

    def __repr__(self):
        return self.__str__()
