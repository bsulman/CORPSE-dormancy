
import CORPSE
import numpy
from pylab import *

from dormancy_params import params,params_nodormant

# In case you want to run without cluttering things up with all the csv files
save_csv=False

############### initial biomass values based on our observations #######################
# With a bulk density of 1.4 g/cm3 and a depth of 10cm, 50ugC/g soil (aprox. TMB observed) is aprox. 7 gC/m2
# if 97% of TMB is dormant before the rewetting (based on our observations)
# initial dormant biomass = 0.97*7 = 6.79. According to CORPSE biomass is in KgC/m2 (although is ploted in gC/m2), so initial value: 0.00679
# initial active biomass = 0.03*7 = 0.21. According to CORPSE biomass is in KgC/m2 (although is ploted in gC/m2), so initial value: 0.00021

# BNS:
# Soils are 1-6% SOM. Assuming 50% carbon, that's 0.5-3% C. For same bulk density (1.4 g/cm3), and 10cm depth, that is
# 700 - 4200 gC/m2 or 0.7-4.2 kgC/m2


initialFast = 0.011
initialSlow = 1.5
initialNecromass = 0.005

initiallitterC = [initialFast,initialSlow,initialNecromass]

dt=1.0/(365.0*24*2)  # Half-hourly resolution
simlength=(6.0+48.0)/(365.0*24)  # In years
nsteps=int(simlength/dt)

inputs=zeros((nsteps,3))  # Zero inputs to start
inputs_control=inputs.copy()

litterC_control_20=zeros((nsteps,3))
protectedC_control_20=zeros((nsteps,3))
microbeC_control_20=zeros(nsteps)
dmicrobeC_control_20=zeros(nsteps)
CO2_control_20=zeros(nsteps)
CO2prod_control_20=zeros(nsteps)

litterC_nodormant_20=zeros((nsteps,3))
protectedC_nodormant_20=zeros((nsteps,3))
microbeC_nodormant_20=zeros(nsteps)
dmicrobeC_nodormant_20=zeros(nsteps)
CO2_nodormant_20=zeros(nsteps)
CO2prod_nodormant_20=zeros(nsteps)

litterC_dormant_20=zeros((nsteps,3))
protectedC_dormant_20=zeros((nsteps,3))
microbeC_dormant_20=zeros(nsteps)
dmicrobeC_dormant_20=zeros(nsteps)
CO2_dormant_20=zeros(nsteps)
CO2prod_dormant_20=zeros(nsteps)

litterC_control_30=zeros((nsteps,3))
protectedC_control_30=zeros((nsteps,3))
microbeC_control_30=zeros(nsteps)
dmicrobeC_control_30=zeros(nsteps)
CO2_control_30=zeros(nsteps)
CO2prod_control_30=zeros(nsteps)

litterC_nodormant_30=zeros((nsteps,3))
protectedC_nodormant_30=zeros((nsteps,3))
microbeC_nodormant_30=zeros(nsteps)
dmicrobeC_nodormant_30=zeros(nsteps)
CO2_nodormant_30=zeros(nsteps)
CO2prod_nodormant_30=zeros(nsteps)

litterC_dormant_30=zeros((nsteps,3))
protectedC_dormant_30=zeros((nsteps,3))
microbeC_dormant_30=zeros(nsteps)
dmicrobeC_dormant_30=zeros(nsteps)
CO2_dormant_30=zeros(nsteps)
CO2prod_dormant_30=zeros(nsteps)

########################################################################################
# c=CORPSE.soil_carbon_cohort(litterC=initiallitterC,activeMicrobeC=initialAMB,dormantMicrobeC=initialDMB,protectedC=[0.1,0.0,0.47],params=params)
from dormancy_params import c_spunup as c
c.set_params(params)
c_control=c.copy()
c_nodormant=c.copy(); c_nodormant.set_params(params_nodormant)
########################################################################################


c_20=c.copy()
c_30=c.copy()
c_nodormant_20=c_nodormant.copy()
c_nodormant_30=c_nodormant.copy()
c_control_20=c_control.copy()
c_control_30=c_control.copy()



starttime=48.0/(365*24)
starttime2=53.0/(365*24)
start=int(starttime/dt)
start2=int(starttime2/dt)

theta=zeros(nsteps)+0.1
theta_control=theta.copy()
theta[start:start2]=0.6  # Paper says wetted to 60% WHC


for step in range(nsteps):

	output_control=c_control_20.update(293.15,theta_control[step],dt) 
	c_control_20.check_validity()

	c_control_20.add_carbon(inputs_control[step,:]*dt)

	litterC_control_20[step,:]=c_control_20.litterC
	protectedC_control_20[step,:]=c_control_20.protectedC
	microbeC_control_20[step]=c_control_20.activeMicrobeC
	dmicrobeC_control_20[step]=c_control_20.dormantMicrobeC
	CO2_control_20[step]=c_control_20.CO2
	CO2prod_control_20[step]=output_control['CO2prod']/dt
	output=c_20.update(293.15,theta[step],dt)                          
	c_20.check_validity()
	c_20.add_carbon(inputs[step,:]*dt)

	litterC_dormant_20[step,:]=c_20.litterC
	protectedC_dormant_20[step,:]=c_20.protectedC
	microbeC_dormant_20[step]=c_20.activeMicrobeC
	dmicrobeC_dormant_20[step]=c_20.dormantMicrobeC
	CO2_dormant_20[step]=c_20.CO2
	CO2prod_dormant_20[step]=output['CO2prod']/dt

	output_nodormant=c_nodormant_20.update(293.15,theta[step],dt)       
	c_nodormant_20.check_validity()

	c_nodormant_20.add_carbon(inputs[step,:]*dt)

	litterC_nodormant_20[step,:]=c_nodormant_20.litterC
	protectedC_nodormant_20[step,:]=c_nodormant_20.protectedC
	microbeC_nodormant_20[step]=c_nodormant_20.activeMicrobeC
	dmicrobeC_nodormant_20[step]=c_nodormant_20.dormantMicrobeC
	CO2_nodormant_20[step]=c_nodormant_20.CO2
	CO2prod_nodormant_20[step]=output_nodormant['CO2prod']/dt


	output_control=c_control_30.update(303.15,theta_control[step],dt)
	c_control_30.check_validity()

	c_control_30.add_carbon(inputs_control[step,:]*dt)

	litterC_control_30[step,:]=c_control_30.litterC
	protectedC_control_30[step,:]=c_control_30.protectedC
	microbeC_control_30[step]=c_control_30.activeMicrobeC
	dmicrobeC_control_30[step]=c_control_30.dormantMicrobeC
	CO2_control_30[step]=c_control_30.CO2
	CO2prod_control_30[step]=output_control['CO2prod']/dt
	output=c_30.update(303.15,theta[step],dt)
	c_30.check_validity()
	c_30.add_carbon(inputs[step,:]*dt)

	litterC_dormant_30[step,:]=c_30.litterC
	protectedC_dormant_30[step,:]=c_30.protectedC
	microbeC_dormant_30[step]=c_30.activeMicrobeC
	dmicrobeC_dormant_30[step]=c_30.dormantMicrobeC
	CO2_dormant_30[step]=c_30.CO2
	CO2prod_dormant_30[step]=output['CO2prod']/dt
	output_nodormant=c_nodormant_30.update(303.15,theta[step],dt)
	c_nodormant_30.check_validity()

	c_nodormant_30.add_carbon(inputs[step,:]*dt)

	litterC_nodormant_30[step,:]=c_nodormant_30.litterC
	protectedC_nodormant_30[step,:]=c_nodormant_30.protectedC
	microbeC_nodormant_30[step]=c_nodormant_30.activeMicrobeC
	dmicrobeC_nodormant_30[step]=c_nodormant_30.dormantMicrobeC
	CO2_nodormant_30[step]=c_nodormant_30.CO2
	CO2prod_nodormant_30[step]=output_nodormant['CO2prod']/dt



# Create CVS files
if save_csv:
#	numpy.savetxt("t.csv",t, delimiter=",")
#	numpy.savetxt("WHC.csv",theta, delimiter=",")
	numpy.savetxt("CO2_dor_20.csv",((CO2prod_dormant_20*1e6/(24*365))/1000), delimiter=",") # in g C m-2 h-1
	numpy.savetxt("CO2_nodor_20.csv",((CO2prod_nodormant_20*1e6/(24*365))/1000), delimiter=",") # in g C m-2 h-1
	numpy.savetxt("AMB_20.csv",microbeC_dormant_20*1000, delimiter=",")
	numpy.savetxt("log10AMB_20.csv",log10((microbeC_dormant_20*1000)+1), delimiter=",")
	numpy.savetxt("DMB_20.csv",dmicrobeC_dormant_20*1000, delimiter=",")
	numpy.savetxt("TMB_dor_20.csv",(microbeC_dormant_20*1000+dmicrobeC_dormant_20*1000), delimiter=",")
	numpy.savetxt("TMB_nodor_20.csv",(microbeC_nodormant_20*1000+dmicrobeC_nodormant_20*1000), delimiter=",")
	numpy.savetxt("FAMB_20.csv",(microbeC_dormant_20/(microbeC_dormant_20+dmicrobeC_dormant_20)), delimiter=",")

	numpy.savetxt("CO2_dor_30.csv",((CO2prod_dormant_30*1e6/(24*365))/1000), delimiter=",") # in g C m-2 h-1
	numpy.savetxt("CO2_nodor_30.csv",((CO2prod_nodormant_30*1e6/(24*365))/1000), delimiter=",") # in g C m-2 h-1
	numpy.savetxt("AMB_30.csv",microbeC_dormant_30*1000, delimiter=",")
	numpy.savetxt("log10AMB_30.csv",log10((microbeC_dormant_30*1000)+1), delimiter=",")
	numpy.savetxt("DMB_30.csv",dmicrobeC_dormant_30*1000, delimiter=",")
	numpy.savetxt("TMB_dor_30.csv",(microbeC_dormant_30*1000+dmicrobeC_dormant_30*1000), delimiter=",")
	numpy.savetxt("TMB_nodor_30.csv",(microbeC_nodormant_30*1000+dmicrobeC_nodormant_30*1000), delimiter=",")
	numpy.savetxt("FAMB_30.csv",(microbeC_dormant_30/(microbeC_dormant_30+dmicrobeC_dormant_30)), delimiter=",")


# With bulk density of 1.4 g/cm3, conversion of kgC/m3 is 1/(1.4e-3kg/100e-3m3)
m_g_conversion=1.0/(1.4e-3/(100**2))

plotstart=int(47/(dt*24*365))
t=arange(nsteps-plotstart)*dt*(24*365)

figure(1);clf()
subplot(311)
plot(t,CO2prod_dormant_20[plotstart:]*1e6/(24*365*1000),'b-',label='Dormancy')
plot(t,CO2prod_control_20[plotstart:]*1e6/(24*365*1000),'b--',label='Control')
plot(t,CO2prod_nodormant_20[plotstart:]*1e6/(24*365*1000),'b:',label='No dormancy')	#AS
plot(t,CO2prod_dormant_30[plotstart:]*1e6/(24*365*1000),'r-')
plot(t,CO2prod_control_30[plotstart:]*1e6/(24*365*1000),'r--')
plot(t,CO2prod_nodormant_30[plotstart:]*1e6/(24*365*1000),'r:')	#AS
xlabel('Time (hours)')
ylabel('CO$_2$ production (gC m$^{-2}$ hour$^{-1}$)')
legend(loc='best',fontsize='medium').get_frame().set_alpha(0.5)
ylim(0,0.5)


# figure(2);clf()
subplot(312)

plot(t,microbeC_dormant_20[plotstart:]*1000,'b-',label='Active Mic')	  # AS changed active from blue 'b' to red 'r'
plot(t,dmicrobeC_dormant_20[plotstart:]*1000,'b--',label='Dormant Mic')	# AS changed dormant from red 'r' to blue 'b'
# plot(t,protectedC_dormant.sum(axis=1),'m-',label='Protected')
plot(t,microbeC_dormant_20[plotstart:]*1000+dmicrobeC_dormant_20[plotstart:]*1000,'b:',label='Total Mic')

plot(t,microbeC_dormant_30[plotstart:]*1000,'r-')	  # AS changed active from blue 'b' to red 'r'
plot(t,dmicrobeC_dormant_30[plotstart:]*1000,'r--')	# AS changed dormant from red 'r' to blue 'b'
# plot(t,protectedC_dormant.sum(axis=1),'m-',label='Protected')
plot(t,microbeC_dormant_30[plotstart:]*1000+dmicrobeC_dormant_30[plotstart:]*1000,'r:')

# plot(t,microbeC_control*1000,'r--')
# plot(t,dmicrobeC_control*1000,'b--')
# plot(t,microbeC_control*1000+dmicrobeC_control*1000,'k--')

plot(t,microbeC_nodormant_20[plotstart:]*1000,'b-.',label='No dormancy')							#AS
plot(t,microbeC_nodormant_20[plotstart:]*1000+dmicrobeC_nodormant_20[plotstart:]*1000,'b-.')   #AS
plot(t,microbeC_nodormant_30[plotstart:]*1000,'r-.')							#AS
plot(t,microbeC_nodormant_30[plotstart:]*1000+dmicrobeC_nodormant_30[plotstart:]*1000,'r-.')   #AS

xlabel('Time (hours)')
ylabel('Carbon pools (gC m$^{-2}$)')

legend(loc='best',fontsize='medium').get_frame().set_alpha(0.5)

subplot(313)
plot(t,microbeC_dormant_20[plotstart:]/(microbeC_dormant_20[plotstart:]+dmicrobeC_dormant_20[plotstart:]),'b-')
plot(t,microbeC_nodormant_20[plotstart:]/(microbeC_nodormant_20[plotstart:]+dmicrobeC_nodormant_20[plotstart:]),'b:')
plot(t,microbeC_control_20[plotstart:]/(microbeC_control_20[plotstart:]+dmicrobeC_control_20[plotstart:]),'b--')

plot(t,microbeC_dormant_30[plotstart:]/(microbeC_dormant_30[plotstart:]+dmicrobeC_dormant_30[plotstart:]),'r-')
plot(t,microbeC_nodormant_30[plotstart:]/(microbeC_nodormant_30[plotstart:]+dmicrobeC_nodormant_30[plotstart:]),'r:')
plot(t,microbeC_control_30[plotstart:]/(microbeC_control_30[plotstart:]+dmicrobeC_control_30[plotstart:]),'r--')
ylim(0,0.2)

xlabel('Time (hours)')
ylabel('Active fraction')

draw()
#
# figure(3);clf()
# plot(t,litterC_dormant[:,0],'b-',label='Fast')
# plot(t,litterC_control[:,0],'b--')
# plot(t,litterC_nodormant[:,0],'b:')
# plot(t,litterC_dormant[:,1],'g-',label='Slow')
# plot(t,litterC_control[:,1],'g--')
# plot(t,litterC_nodormant[:,1],'g:')
# plot(t,litterC_dormant[:,2],'r-',label='Necromass')
# plot(t,litterC_control[:,2],'r--')
# plot(t,litterC_nodormant[:,2],'r:')
#
# title('Substrate C pools')
# xlabel('Time (hours)')
# ylabel('Carbon pools (kgC m$^{-2}$)')
# legend(loc='best')
#
# draw()

show()
