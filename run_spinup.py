
import CORPSE
import numpy
from pylab import *

from dormancy_params import params,params_nodormant

save_csv=False

dt=1.0/365 # Daily time step
simlength=30.0 # years
nsteps=int(simlength/dt)

maxT=25.0 # These are based on normal annual pattern of daily air temperature at MMSF
minT=-10.0
# T=zeros(nsteps)+273.15+10
T=(-cos(arange(nsteps)*2*pi/365.0)+1)*(maxT-minT)/2.0+minT+273.15   # Sinusoid curve with annual time scale
theta=zeros(nsteps)+0.6
# Make inputs equivalent to 120 mgC/m2/hour at constant rate
# So we can calculate things at equilibrium for those conditions
inputs_total=120e-6*24*365*0.5
inputs=zeros((nsteps,3))
inputs[:,0]=0.1*inputs_total;inputs[:,1]=0.9*inputs_total;inputs[:,2]=0.0*inputs_total

# Do last four days at 20C with some wetting just to check how things look
inputs[-4:,:]=0.0
T[-4:]=293.15
theta[-4:-2]=0.1
theta[-2:]=0.6

initialAMB = 0.00021
initialDMB = 0.00679
initialTMB = initialAMB + initialDMB

initialFast = 0.011
initialSlow = 1.5
initialNecromass = 0.005

initiallitterC = [initialFast,initialSlow,initialNecromass]
initialProtectedC = [0.1,0.0,0.47]

c=CORPSE.soil_carbon_cohort(litterC=initiallitterC,activeMicrobeC=initialAMB,dormantMicrobeC=initialDMB,protectedC=initialProtectedC,params=params)

c_nodormant=CORPSE.soil_carbon_cohort(litterC=initiallitterC,activeMicrobeC=initialTMB,dormantMicrobeC=0.0,protectedC=initialProtectedC,params=params_nodormant)

litterC_nodormant=zeros((nsteps,3))
protectedC_nodormant=zeros((nsteps,3))
microbeC_nodormant=zeros(nsteps)
dmicrobeC_nodormant=zeros(nsteps)
CO2_nodormant=zeros(nsteps)
CO2prod_nodormant=zeros(nsteps)

litterC_dormant=zeros((nsteps,3))
protectedC_dormant=zeros((nsteps,3))
microbeC_dormant=zeros(nsteps)
dmicrobeC_dormant=zeros(nsteps)
CO2_dormant=zeros(nsteps)
CO2prod_dormant=zeros(nsteps)

for step in range(nsteps):

	output=c.update(T[step],theta[step],dt)
	c.check_validity()
	c.add_carbon(inputs[step,:]*dt)

	litterC_dormant[step,:]=c.litterC
	protectedC_dormant[step,:]=c.protectedC
	microbeC_dormant[step]=c.activeMicrobeC
	dmicrobeC_dormant[step]=c.dormantMicrobeC
	CO2_dormant[step]=c.CO2
	CO2prod_dormant[step]=output['CO2prod']/dt


	output_nodormant=c_nodormant.update(T[step],theta[step],dt)
	c_nodormant.check_validity()
	c_nodormant.add_carbon(inputs[step,:]*dt)

	litterC_nodormant[step,:]=c_nodormant.litterC
	protectedC_nodormant[step,:]=c_nodormant.protectedC
	microbeC_nodormant[step]=c_nodormant.activeMicrobeC
	dmicrobeC_nodormant[step]=c_nodormant.dormantMicrobeC
	CO2_nodormant[step]=c_nodormant.CO2
	CO2prod_nodormant[step]=output_nodormant['CO2prod']/dt


# Create CVS files
if save_csv:
	numpy.savetxt("t.csv",t, delimiter=",")
	numpy.savetxt("WHC_spinup.csv",theta, delimiter=",")
	numpy.savetxt("CO2_dor_spinup.csv",(CO2prod_dormant*1e6/(24*365)), delimiter=",")
	numpy.savetxt("CO2_nodor_spinup.csv",(CO2prod_nodormant*1e6/(24*365)), delimiter=",")
	numpy.savetxt("AMB_spinup.csv",microbeC_dormant*1000, delimiter=",")
	numpy.savetxt("log10AMB_spinup.csv",log10((microbeC_dormant*1000)+1), delimiter=",")
	numpy.savetxt("DMB_spinup.csv",dmicrobeC_dormant*1000, delimiter=",")
	numpy.savetxt("TMB_dor_spinup.csv",(microbeC_dormant*1000+dmicrobeC_dormant*1000), delimiter=",")
	numpy.savetxt("TMB_nodor_spinup.csv",(microbeC_nodormant*1000+dmicrobeC_nodormant*1000), delimiter=",")
	numpy.savetxt("FAMB_spinup.csv",(microbeC_dormant/(microbeC_dormant+dmicrobeC_dormant)), delimiter=",")

t=arange(nsteps)*dt

m_g_conversion=1.0/(1.4e-3/(100**2))

figure(1);clf()
subplot(311)
plot(t,CO2prod_dormant*1e6/(24*365),'k-',label='Dormancy')
plot(t,CO2prod_nodormant*1e6/(24*365),'k:',label='No dormancy')	#AS
xlabel('Time (hours)')
ylabel('CO$_2$ production (mgC m$^{-2}$ hour$^{-1}$)')
legend(loc='best')

subplot(312)
# plot(t,litterC_dormant[:,0],'b-',label='Fast')
# plot(t,litterC_dormant[:,1],'g-',label='Slow')
# plot(t,litterC_dormant[:,2],'r-',label='Dead mic')

plot(t,microbeC_dormant*1000,'r-',label='Active Mic')	  # AS changed active from blue 'b' to red 'r'
plot(t,dmicrobeC_dormant*1000,'b-',label='Dormant Mic')	# AS changed dormant from red 'r' to blue 'b'
# plot(t,protectedC_dormant.sum(axis=1),'m-',label='Protected')
plot(t,microbeC_dormant*1000+dmicrobeC_dormant*1000,'k-',label='Total Mic')

# plot(t,litterC_nodormant[:,0],'b:')
# plot(t,litterC_nodormant[:,1],'g:')
# plot(t,litterC_nodormant[:,2],'r:')
plot(t,microbeC_nodormant*1000,'b:')							#AS
# plot(t,dmicrobeC_nodormant*1000,'k:')
# plot(t,protectedC_nodormant.sum(axis=1),'m:')
plot(t,microbeC_nodormant*1000+dmicrobeC_nodormant*1000,'k:')   #AS

xlabel('Time (years)')
ylabel('Carbon pools (gC m$^{-2}$)')

legend(loc='best').get_frame().set_alpha(0.5)

subplot(313)
plot(t,microbeC_dormant/(microbeC_dormant+dmicrobeC_dormant),'k-')
plot(t,microbeC_nodormant/(microbeC_nodormant+dmicrobeC_nodormant),'k:')
ylim(0,1.0)

xlabel('Time (years)')
ylabel('Active fraction')

draw()

figure(2);clf()
plot(t,litterC_dormant[:,0],'b-',label='Fast')
plot(t,litterC_nodormant[:,0],'b:')
plot(t,litterC_dormant[:,1],'g-',label='Slow')
plot(t,litterC_nodormant[:,1],'g:')
plot(t,litterC_dormant[:,2],'r-',label='Necromass')
plot(t,litterC_nodormant[:,2],'r:')

title('Substrate C pools')
xlabel('Time (years)')
ylabel('Carbon pools (kgC m$^{-2}$)')
legend(loc='best')

draw()


# Run for two more days at 20C or 30C
dt=1.0/(24.0*365.0) # Hourly resolution
nsteps=48+24

c_20=c.copy()
c_30=c.copy()
c_nodormant_20=c_nodormant.copy()
c_nodormant_30=c_nodormant.copy()

CO2prod_20=zeros(nsteps)
CO2prod_30=zeros(nsteps)
CO2prod_nodormant_20=zeros(nsteps)
CO2prod_nodormant_30=zeros(nsteps)

microbeC_20=zeros(nsteps)
dmicrobeC_20=zeros(nsteps)
microbeC_30=zeros(nsteps)
dmicrobeC_30=zeros(nsteps)
microbeC_nodormant_20=zeros(nsteps)
dmicrobeC_nodormant_20=zeros(nsteps)
microbeC_nodormant_30=zeros(nsteps)
dmicrobeC_nodormant_30=zeros(nsteps)

theta_lab=zeros(nsteps)
theta_lab[:48]=0.1
theta_lab[48:]=0.6

for step in range(nsteps):

	output=c_20.update(20+273.15,theta_lab[step],dt)
	c_20.check_validity()
	CO2prod_20[step]=output['CO2prod']/dt
	microbeC_20[step]=c_20.activeMicrobeC
	dmicrobeC_20[step]=c_20.dormantMicrobeC

	output=c_30.update(30+273.15,theta_lab[step],dt)
	c_30.check_validity()
	CO2prod_30[step]=output['CO2prod']/dt
	microbeC_30[step]=c_30.activeMicrobeC
	dmicrobeC_30[step]=c_30.dormantMicrobeC


	output=c_nodormant_20.update(20+273.15,theta_lab[step],dt)
	c_nodormant_20.check_validity()
	CO2prod_nodormant_20[step]=output['CO2prod']/dt
	microbeC_nodormant_20[step]=c_nodormant_20.activeMicrobeC
	dmicrobeC_nodormant_20[step]=c_nodormant_20.dormantMicrobeC

	output=c_nodormant_30.update(30+273.15,theta_lab[step],dt)
	c_nodormant_30.check_validity()
	CO2prod_nodormant_30[step]=output['CO2prod']/dt
	microbeC_nodormant_30[step]=c_nodormant_30.activeMicrobeC
	dmicrobeC_nodormant_30[step]=c_nodormant_30.dormantMicrobeC

figure(3);clf()
t2=arange(nsteps)
subplot(311)
plot(t2,CO2prod_20*1e6/(24*365),'b-',label='Dormancy (20C)')
plot(t2,CO2prod_30*1e6/(24*365),'r-',label='Dormancy (30C)')
plot(t2,CO2prod_nodormant_20*1e6/(24*365),'b--',label='No Dormancy (20C)')
plot(t2,CO2prod_nodormant_30*1e6/(24*365),'r--',label='No Dormancy (30C)')
xlabel('Time (hours)')
ylabel('CO$_2$ production (mgC m$^{-2}$ hour$^{-1}$)')
legend(loc='best')
title('Lab conditions')

subplot(312)

plot(t2,log10(microbeC_20*1000),'b-',label='Active Mic')
plot(t2,log10(dmicrobeC_20*1000),'b--',label='Dormant Mic')
plot(t2,log10(microbeC_20*1000+dmicrobeC_20*1000),'b:',label='Total Mic')

plot(t2,log10(microbeC_30*1000),'r-')
plot(t2,log10(dmicrobeC_30*1000),'r--')
plot(t2,log10(microbeC_30*1000+dmicrobeC_30*1000),'r:')

plot(t2,log10(microbeC_nodormant_20*1000),'b-.')
plot(t2,log10(microbeC_nodormant_30*1000),'r-.')

xlabel('Time (hours)')
ylabel('log10 mic biomass (gC m$^{-2}$)')

legend(loc='best').get_frame().set_alpha(0.5)

subplot(313)
def logit(x):
	return log(x/(1.0-x))

plot(t2,logit(microbeC_20/(microbeC_20+dmicrobeC_20)),'b-')
plot(t2,logit(microbeC_30/(microbeC_30+dmicrobeC_30)),'r-')
plot(t2,logit(microbeC_nodormant_20/(microbeC_nodormant_20+dmicrobeC_nodormant_20)),'b--')
plot(t2,logit(microbeC_nodormant_30/(microbeC_nodormant_30+dmicrobeC_nodormant_30)),'b--')
# ylim(0,1.0)

xlabel('Time (hours)')
ylabel('logit(Active fraction)')

draw()

show()
