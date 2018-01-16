import nest
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

# display recordables for illustration
print 'iaf_cond_alpha recordables: ', nest.GetDefaults('iaf_cond_alpha')['recordables']
# create neuron and multimeter
n = nest.Create('iaf_cond_alpha',  params = {'tau_syn_ex': 1.0, 'V_reset': -70.0})
m = nest.Create('multimeter', params = {'withtime': True, 'interval': 0.1, 'record_from': ['V_m', 'g_ex', 'g_in']})
s = nest.Create('spike_detector', params = {'withtime': True, 'withgid': True })
s_input_ex = nest.Create('spike_detector', params = {'withtime': True, 'withgid': True })
s_input_inh = nest.Create('spike_detector', params = {'withtime': True, 'withgid': True })
# Create spike generators and connect

stimes=[10.]
stimes+=[20.]
stimes+=[30.]
stimes+=[i*2.0+40. for i in range(30)]
for i in range(2):
    stimes+=[i*2.0+140.]
for i in range(3):
    stimes+=[i*2.0+150.]
for i in range(4):
    stimes+=[i*2.0+160.]
for i in range(5):
    stimes+=[i*2.0+220.]
#for i in range(5):
#    stimes+=[i*2.0+260.}
#stimes+=[261.]
stimes+=[260.]
gex = nest.Create('spike_generator', params = {'spike_times': np.array(stimes)}) #np.array([10.0, 12.0, 14.0, 16.0, 18.0, 26.0, 27.0, 28.0, 29.0, 35., 50.0, 52.0, 54.0, 56.0, 58.0, 66.0, 67.0, 68.0, 69.0, 70.0, 72.0, 79.0, 88.0])})
    #np.array([i*2.0+10. for i in range(35)])})
    
#stimes=[31.]
stimes=[i*10.0+40. for i in range(7)]
stimes+=[i*10.0+150. for i in range(4)]
#for i in range(3):
#    stimes+=[i*10.0+180.]
stimes+=[232.]
stimes+=[260.]

gin = nest.Create('spike_generator',  params = {'spike_times': np.array(stimes)}) #np.array([10.0, 25.0, 40.0])})
nest.Connect(gex, n, 'one_to_one', {'weight':  50.0, 'delay' : 12.}) # excitatory #40.
nest.Connect(gin, n, 'one_to_one', {'weight': -400.0}) # inhibitory
nest.Connect(gex, s_input_ex)
nest.Connect(gin, s_input_inh)
nest.Connect(m, n)
nest.Connect(n, s)
# simulate
nest.Simulate(300)
# obtain and display data
events = nest.GetStatus(m)[0]['events']
t = events['times'];

f, (a0, a1) = plt.subplots(2, 1, sharex='col', gridspec_kw = {'height_ratios':[5, 1]}, figsize=(6.5,2.5))

a0.plot(t, (events['V_m']+70.)/15., 'k')
a0.plot([0,300], [1.,1.], 'k--')
a0.set_xlim([0, 100])
a0.set_ylim([-0.1, 1.1])
a0.set_ylabel('Membrane potential\n [arbitray units]', fontsize=8)
a0.tick_params(labelsize=8)

dSD = nest.GetStatus(s,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
a1.plot(ts, evs, "k.")

dSD = nest.GetStatus(s_input_ex,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
a1.plot(ts, evs/10.+0.3, ".", c="gray")

dSD = nest.GetStatus(s_input_inh,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
a1.plot(ts, evs/10.+0.25, ".", c="red")

#pl.plot(t, events['g_ex'], t, events['g_in'])
#a1.set_xlim([0, 60])
#a1.set_ylim([0.9, 1.1])
#a1.set_xlabel('Time [ms]')
#a1.set_ylabel('Input')
#a1.set_yticks([])

#pl.plot(t, events['g_ex'], t, events['g_in'])
a1.set_xlim([0, 300])
a1.set_ylim([0.85, 1.05])
a1.set_xlabel('Time [ms]', fontsize=8)
a1.set_ylabel('Spikes', fontsize=8)
a1.set_yticks([])
a1.tick_params(labelsize=8)
#pl.ylabel('Synaptic conductance [nS]')
plt.legend(('Output', 'Input ex.', 'Input inh.'), loc='upper center', bbox_to_anchor=(0.5,1.9), frameon=False, fancybox=False, ncol=3, fontsize=8)
plt.tight_layout()
plt.savefig('img/ana/all_mo_det_multi_large.pgf')
plt.savefig('img/ana/all_mo_det_multi_large.pdf')
plt.show()