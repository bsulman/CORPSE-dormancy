import CORPSE
import numpy
from pylab import *

from dormancy_params import params,params_nodormant

save_csv=False

num_cycles=15
cycle_length=50/(365.0) #  20/(365.0): Between wetting events (currently 10 days, in years units)

dt=1.0/(365.0*24)  # Hourly resolution
simlength=cycle_length*num_cycles  # In years
nsteps=int(simlength/dt)

inputs=zeros((nsteps,3))  # Zero inputs to start
inputs_control=inputs.copy()
theta=zeros(nsteps)+0.1
T=273.15+20

theta_decay_time = 5.0/(365.0) # Exponential decay rate

iter_time=((arange(nsteps)-cycle_length*0.25/dt)%(cycle_length/dt))*dt
theta=0.1+0.5*exp(-iter_time/theta_decay_time)
theta_control=zeros(nsteps)+0.1

litterC_control=zeros((nsteps,3))
protectedC_control=zeros((nsteps,3))
microbeC_control=zeros(nsteps)
dmicrobeC_control=zeros(nsteps)
CO2_control=zeros(nsteps)
CO2prod_control=zeros(nsteps)

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


from dormancy_params import c_spunup as c
c.set_params(params)
c_control=c.copy()
c_nodormant=c.copy()
c_nodormant.set_params(params_nodormant)


for step in range(nsteps):

    output_control=c_control.update(T,theta_control[step],dt)
    c_control.check_validity()

    c_control.add_carbon(inputs_control[step,:]*dt)

    litterC_control[step,:]=c_control.litterC
    protectedC_control[step,:]=c_control.protectedC
    microbeC_control[step]=c_control.activeMicrobeC
    dmicrobeC_control[step]=c_control.dormantMicrobeC
    CO2_control[step]=c_control.CO2
    CO2prod_control[step]=output_control['CO2prod']/dt


    output=c.update(T,theta[step],dt)
    c.check_validity()
    c.add_carbon(inputs[step,:]*dt)

    litterC_dormant[step,:]=c.litterC
    protectedC_dormant[step,:]=c.protectedC
    microbeC_dormant[step]=c.activeMicrobeC
    dmicrobeC_dormant[step]=c.dormantMicrobeC
    CO2_dormant[step]=c.CO2
    CO2prod_dormant[step]=output['CO2prod']/dt

    output_nodormant=c_nodormant.update(T,theta[step],dt)
    c_nodormant.check_validity()

    c_nodormant.add_carbon(inputs[step,:]*dt)

    litterC_nodormant[step,:]=c_nodormant.litterC
    protectedC_nodormant[step,:]=c_nodormant.protectedC
    microbeC_nodormant[step]=c_nodormant.activeMicrobeC
    dmicrobeC_nodormant[step]=c_nodormant.dormantMicrobeC
    CO2_nodormant[step]=c_nodormant.CO2
    CO2prod_nodormant[step]=output_nodormant['CO2prod']/dt


t=arange(nsteps)*dt


# Create CVS files

if save_csv:
#	numpy.savetxt("t.csv",t, delimiter=",")
#	numpy.savetxt("WHC_WC.csv",theta, delimiter=",")
#	CO2_dor=numpy.savetxt("CO2_dor_WC.csv",(CO2prod_dormant*1e6/(24*365*1000)), delimiter=",") # in g C m-2 h-1
#	CO2_nodor=numpy.savetxt("CO2_nodor_WC.csv",(CO2prod_nodormant*1e6/(24*365*1000)), delimiter=",") # in g C m-2 h-1
#	numpy.savetxt("AMB_WC.csv",microbeC_dormant*1000, delimiter=",")
#	numpy.savetxt("log10AMB_WC.csv",log10((microbeC_dormant*1000)+1), delimiter=",")
#	numpy.savetxt("DMB_WC.csv",dmicrobeC_dormant*1000, delimiter=",")
#	numpy.savetxt("TMB_dor_WC.csv",(microbeC_dormant*1000+dmicrobeC_dormant*1000), delimiter=",")
#	numpy.savetxt("TMB_nodor_WC.csv",(microbeC_nodormant*1000+dmicrobeC_nodormant*1000), delimiter=",")
#	numpy.savetxt("FAMB_WC.csv",(microbeC_dormant/(microbeC_dormant+dmicrobeC_dormant)), delimiter=",")
#	numpy.savetxt("litterC_slow_dor.csv",litterC_dormant[:,1], delimiter=",")
#	numpy.savetxt("litterC_slow_nodor.csv",litterC_nodormant[:,1], delimiter=",")
	numpy.savetxt("table.csv", np.transpose([(CO2prod_dormant*1e6/(24*365*1000)), # CO2 dormancy model
                                          (CO2prod_nodormant*1e6/(24*365*1000)),   # CO2 no-dormancy model
                                          (microbeC_dormant*1000+dmicrobeC_dormant*1000),  # TMB dormancy model
                                          (microbeC_nodormant*1000+dmicrobeC_nodormant*1000),  # TMB no-dormancy model
                                          (microbeC_dormant/(microbeC_dormant+dmicrobeC_dormant)),  # FAMB
                                          (litterC_dormant[:,1]),   # Slow substrate C pool dormancy model
                                          (litterC_nodormant[:,1])]), delimiter=",") # Slow substrate C pool no-dormancy model


figure(1);clf()
subplot(311)
plot(t,CO2prod_dormant*1e6/(24*365*1000),'k-',label='Dormancy')
plot(t,CO2prod_nodormant*1e6/(24*365*1000),'k:',label='No dormancy')	#AS
xlabel('Time (hours)')
ylabel('CO$_2$ production (gC m$^{-2}$ hour$^{-1}$)')
legend(loc='best',fontsize='small').get_frame().set_alpha(0.6)

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
ylim(0,0.05)

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
show()
