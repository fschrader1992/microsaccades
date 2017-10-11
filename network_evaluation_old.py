#-------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------SAVING-AND-PRINTING-THE-OUTPUT
now = datetime.datetime.now()

store_sp = []
store_sp_200=[]
s_gids = []

fig = plt.figure()
fig.set_size_inches(6,5)
ax = plt.subplot(111)

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/midgets_'+str(handle_name)+'.txt','w+')

firstq = True
#poisson midget cells
h_rates=[]
#for q in range(1,121):
for pos in gm_pos:
    s = tp.FindNearestElement(out_m,[pos[0],pos[1]]) #[7.5,float(q/2.)]) #[float(q/2.),7.5]
    #s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200+=[h200]
        #print q,s,sp
        #qa = [q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'r.',label='spikes')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'r.',label='_nolegend_')
        
        if h>1:
            h_rate_loc=0
            for u in range(h-1):
                 h_rate_loc+=sp[u+1]-sp[u]
            h_rate_loc /= (h-1)
            h_rates +=[h_rate_loc]
        
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

m_in_max_spikes = max(store_sp)
if m_in_max_spikes == 0:
    m_in_max_spikes = [0]
print m_in_max_spikes

m_in_max_spikes_200 = [max(store_sp_200)]
if m_in_max_spikes_200 == 0:
    m_in_max_spikes_200 = [0]
print str(m_in_max_spikes_200)+' HERE'

as_file = open('data/mo_det_cal/average_spiking_I_e_f0.txt','a+')
as_file.write(str(I_E)+'\t'+str(STDI)+'\t'+str(np.mean(h_rates))+'\t'+str(np.std(h_rates))+'\n')
as_file.close()

all_movement_spikes = 0.
all_movement_spikes_200 = 0.

#left-right-part----------------------------------------------------------------------------------------------------------------

#pyl.figure()
#V_E = nest.GetStatus(s_gids, 'rate')
#pyl.hist(V_E, bins = 100)
#pyl.show()

#midget motion detectors complete
#ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')

#midget motion detectors left
store_sp = []
store_sp_200=[]
s_gids = []

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/midgets_md_right_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,121):
for pos in gm_pos:
    s = tp.FindNearestElement(out_m_r_0_right,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        #qa = [q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='leftward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')

    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_right = max(store_sp)
if max_spikes_right == 0:
    max_spikes_right = [0]
print max_spikes_right

max_spikes_right_200 = [max(store_sp_200)]
if max_spikes_right_200 == 0:
    max_spikes_right_200 = [0]


#midget motion detectors left
store_sp = []
store_sp_200=[]
s_gids = []

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/midgets_md_left_'+str(handle_name)+'.txt','w+')

firstq = True
#for pos in gm_pos:
#for q in range(1,121):
for pos in gm_pos:
    s = tp.FindNearestElement(out_m_r_0_left,[pos[0],pos[1]]) #[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        #qa = [q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='rightward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')

    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_left = max(store_sp)
if max_spikes_left == 0:
    max_spikes_left = [0]
print max_spikes_left

max_spikes_left_200 = [max(store_sp_200)]
if max_spikes_left_200 == 0:
    max_spikes_left_200 = [0]


#midget motion detectors complete
#ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')
#store_sp = []
#store_sp_200 = []
#s_gids = []

#firstq = True
#for q in range(1,121):
#    s = tp.FindNearestElement(out_m_r_0,[float(q/2.),7.5])
#    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
#    s_gids += [s]
    
#    dSD = nest.GetStatus(s,keys="events")[0]
#    evs = dSD["senders"]
#    ts = dSD["times"]
#    #pyl.figure(2)
#    #pyl.plot(ts, evs, ".")
#    #pyl.show()
    
#    if evs.any():
#        #ev = nest.GetStatus(mult)[0]['events']
#        #t = ev['times']
#        #r = ev['rate']

#        sp = nest.GetStatus(s)[0]['events']['times']
#        #plt.subplot(221)
#        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
#        store_sp +=[h]
#        h200=0
#        for spike in sp:
#            if spike>200.:
#                h200+=1
#        store_sp_200 +=[h200]
        
#        #qa = [q for i in range(h)]
#        #if firstq == True:
#        #    ax.plot(qa,sp,'g*',label='general motion detectors')
#        #    firstq = False
#        #else:
#        #    ax.plot(qa,sp,'g*',label='_nolegend_')
        
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
#        #plt.plot(t, r, color='b')
#        #plt.step(e[:-1], h , color='b', where='post')
#        #plt.title('PST histogram and firing rates')
#        #plt.ylabel('Spikes per second')

#        #plt.subplot(223)
#        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
#                    #histtype='step', color='b')
#        #plt.title('ISI histogram')
#        #plt.show()
#    else:
#        store_sp +=[0]
#        store_sp_200 +=[0]


#ms_all_file.close()
##plt.show()
#max_spikes_lr = max(store_sp)
#if max_spikes_lr == 0:
#    max_spikes_lr = [0]
#print max_spikes_lr

#max_spikes_lr_200 = [max(store_sp_200)]
#if max_spikes_lr_200 == 0:
#    max_spikes_lr_200 = [0]


#up-down-part--------------------------------------------------------------------------------------------------------------------

#midget motion detectors up
store_sp = []
store_sp_200=[]
s_gids = []

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/midgets_md_up_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,121):
for pos in gm_pos:
    s = tp.FindNearestElement(out_m_r_60_up,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        #qa = [q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='upward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')

    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_up = max(store_sp)
if max_spikes_up == 0:
    max_spikes_up = [0]
print max_spikes_up

max_spikes_up_200 = [max(store_sp_200)]
if max_spikes_up_200 == 0:
    max_spikes_up_200 = [0]


#midget motion detectors down
store_sp = []
store_sp_200=[]
s_gids = []

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/midgets_md_down_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,121):
for pos in gm_pos:
    s = tp.FindNearestElement(out_m_r_60_down,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        #qa = [q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='downward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')

    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_down = max(store_sp)
if max_spikes_down == 0:
    max_spikes_down = [0]
print max_spikes_down

max_spikes_down_200 = [max(store_sp_200)]
if max_spikes_down_200 == 0:
    max_spikes_down_200 = [0]


#midget motion detectors complete
#ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')
#store_sp = []
#store_sp_200 = []
#s_gids = []

#firstq = True
#for q in range(1,121):
#    s = tp.FindNearestElement(out_m_r_60,[7.5,float(q/2.)])
#    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
#    s_gids += [s]
#    
#    dSD = nest.GetStatus(s,keys="events")[0]
#    evs = dSD["senders"]
#    ts = dSD["times"]
#    #pyl.figure(2)
#    #pyl.plot(ts, evs, ".")
#    #pyl.show()
#    
#    if evs.any():
#        #ev = nest.GetStatus(mult)[0]['events']
#        #t = ev['times']
#        #r = ev['rate']#
#
#        sp = nest.GetStatus(s)[0]['events']['times']
#        #plt.subplot(221)
#        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
#        store_sp +=[h]
#        h200=0
#        for spike in sp:
#            if spike>200.:
#                h200+=1
#        store_sp_200 +=[h200]
#        
#        qa = [q for i in range(h)]
#        if firstq == True:
#            ax.plot(qa,sp,'k*',label='general motion detectors')
#            firstq = False
#        else:
#            ax.plot(qa,sp,'k*',label='_nolegend_')
#        
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
#        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
#        #plt.plot(t, r, color='b')
#        #plt.step(e[:-1], h , color='b', where='post')
#        #plt.title('PST histogram and firing rates')
#        #plt.ylabel('Spikes per second')

#        #plt.subplot(223)
#        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
#                    #histtype='step', color='b')
#        #plt.title('ISI histogram')
#        #plt.show()
#    else:
#        store_sp +=[0]
#        store_sp_200 +=[0]
        
        
#max_spikes_ud = max(store_sp)
#if max_spikes_ud == 0:
#    max_spikes_ud = [0]
#print max_spikes_ud

#max_spikes_ud_200 = [max(store_sp_200)]
#if max_spikes_ud_200 == 0:
#    max_spikes_ud_200 = [0]
    
#plot/output-part---------------------------------------------------------------------------------------------------------------
#plt.title('midget spikes for motion detectors')
#plt.xlabel('position of neuron in arcmin')
#plt.ylabel('time $t$ in ms')
#plt.legend()
'''
plt.xlim([0,120])
plt.ylim([0,t_end-t_start])

#ax.set_title('spiking times uniformly distributed with $\sigma_{I_e} = 10$pA')
ax.set_xlabel('position of neuron in arcmin')
#ax.get_xaxis().set_visible(False)
ax.set_ylabel('simulated time $t$ in ms')
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

#ax.legend(loc=4, fancybox=False, fontsize=8, shadow=False)
# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), 
          fancybox=False, fontsize=8, shadow=False, ncol=2)

plt.savefig('img/mo_det_cal/x_t_midgets_'+handle_name+'_90deg.pdf')
plt.savefig('img/mo_det_cal/x_t_midgets_'+handle_name+'_90deg.pgf')
out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '_x_t_midgets'+ str(handle_name) + '.pdf'
plt.savefig(out)
plt.show()
plt.close()
'''

#and maximum spikes for one cell
max_spikes_total = max(max_spikes_left[0],max_spikes_right[0],max_spikes_up[0],max_spikes_down[0])
max_spikes_total_200 = max(max_spikes_left_200[0],max_spikes_right_200[0],max_spikes_up_200[0],max_spikes_down_200[0])

ms_file = open('data/poletti2010/m_max_spikes.txt','a+')
#ms_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(max_spikes[0])+'\t'+str(m_in_max_spikes[0])+'\t'+ str(max_spikes_left[0])+'\t'+ str(max_spikes_right[0])+'\t'+ str(max_spikes_200[0])+'\t'+ str(m_in_max_spikes_200[0])+'\t'+ str(max_spikes_left_200[0])+'\t'+ str(max_spikes_right_200[0])+'\n')
#ms_file = open('data/mo_det_cal/m_max_spikes.txt','a+')
ms_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(m_in_max_spikes[0])+'\t'+ str(max_spikes_left[0])+'\t'+ str(max_spikes_right[0])+'\t'+ str(max_spikes_up[0])+'\t'+ str(max_spikes_down[0])+'\t'+str(max_spikes_total)+'\t'+ str(all_movement_spikes)+'\n')
ms_file.close()

ms_file_200 = open('data/poletti2010/m_max_spikes_200.txt','a+')
#ms_file_200 = open('data/mo_det_cal/m_max_spikes_200.txt','a+')
ms_file_200.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(m_in_max_spikes_200[0])+'\t'+ str(max_spikes_left_200[0])+'\t'+ str(max_spikes_right_200[0])+'\t'+ str(max_spikes_up_200[0])+'\t'+ str(max_spikes_down_200[0])+'\t'+str(max_spikes_total_200)+'\t'+ str(all_movement_spikes_200)+'\n')
ms_file_200.close()

print max_spikes_total


#----------------------------------------------------------------------------------------------------------------------------------
#parasol poisson cells-------------------------------------------------------------------------------------------------------------

store_sp=[]
store_sp_200=[]

fig = plt.figure()
fig.set_size_inches(6,5)
ax = plt.subplot(111)

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p,[pos[0],pos[1]]) #[4.*float(q/2.),7.5])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='spikes')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

#plt.show()
p_in_max_spikes = max(store_sp)
if p_in_max_spikes == 0:
    p_in_max_spikes = [0]
print p_in_max_spikes

p_in_max_spikes_200 = [max(store_sp_200)]
if p_in_max_spikes_200 == 0:
    p_in_max_spikes_200 = [0]
    
    
all_movement_spikes = 0.
all_movement_spikes_200 = 0.

#left-right-part----------------------------------------------------------------------------------------------------------------

#parasolic motion detectors right
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_right_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_0_right,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        qa = [4*q for i in range(h)]
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='leftward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_right = max(store_sp)
if max_spikes_right == 0:
    max_spikes_right = [0]
print max_spikes_right

max_spikes_right_200 = [max(store_sp_200)]
if max_spikes_right_200 == 0:
    max_spikes_right_200 = [0]

#parasolic motion detectors left
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_left_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_0_left,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='rightward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_left = max(store_sp)
if max_spikes_left == 0:
    max_spikes_left = [0]
print max_spikes_left

max_spikes_left_200 = [max(store_sp_200)]
if max_spikes_left_200 == 0:
    max_spikes_left_200 = [0]

#parasolic motion detectors complete
#ps_all_file = open('data/mo_det_cal/p_max_spikes_all.txt','a+')

#store_sp=[]
#store_sp_200=[]

#firstq = True
#for q in range(1,31):
#    s = tp.FindNearestElement(out_p_r_0,[4.*float(q/2.),7.5])
#    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
#    s_gids += [s]
    
#    dSD = nest.GetStatus(s,keys="events")[0]
#    evs = dSD["senders"]
#    ts = dSD["times"]
#    #pyl.figure(2)
#    #pyl.plot(ts, evs, ".")
#    #pyl.show()
    
#    if evs.any():
#        #ev = nest.GetStatus(mult)[0]['events']
#        #t = ev['times']
#        #r = ev['rate']

#        sp = nest.GetStatus(s)[0]['events']['times']
#        #plt.subplot(221)
#        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
#        store_sp +=[h]
#        h200=0
#        for spike in sp:
#            if spike>200.:
#                h200+=1
#        store_sp_200 +=[h200]
#        #qa = [4*q for i in range(h)]
#        #if firstq == True:
#        #    ax.plot(qa,sp,'g*',label='general motion detectors')
#        #    firstq = False
#        #else:
#        #    ax.plot(qa,sp,'g*',label='_nolegend_')
#        #ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
#        #ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
#        #ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
#        #plt.plot(t, r, color='b')
#        #plt.step(e[:-1], h , color='b', where='post')
#        #plt.title('PST histogram and firing rates')
#        #plt.ylabel('Spikes per second')

#        #plt.subplot(223)
#        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
#                    #histtype='step', color='b')
#        #plt.title('ISI histogram')
#        #plt.show()
#    else:
#        store_sp +=[0]
#        store_sp_200 +=[0]

#ps_all_file.close()
##plt.show()
#max_spikes_lr = max(store_sp)
#if max_spikes_lr == 0:
#    max_spikes_lr = [0]
#print max_spikes_lr

#max_spikes_lr_200 = [max(store_sp_200)]
#if max_spikes_lr_200 == 0:
#    max_spikes_lr_200 = [0]
##print max_spikes_200


#up-down-part-60--------------------------------------------------------------------------------------------------------------------

#parasolic motion detectors up 60
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_up_60_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_60_up,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='120 downward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_up_60 = max(store_sp)
if max_spikes_up_60 == 0:
    max_spikes_up_60 = [0]
print max_spikes_up_60

max_spikes_up_60_200 = [max(store_sp_200)]
if max_spikes_up_60_200 == 0:
    max_spikes_up_60_200 = [0]

#parasolic motion detectors down 60
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_down_60_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_60_down,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='120 downward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_down_60 = max(store_sp)
if max_spikes_down_60 == 0:
    max_spikes_down_60 = [0]
print max_spikes_down_60

max_spikes_down_60_200 = [max(store_sp_200)]
if max_spikes_down_60_200 == 0:
    max_spikes_down_60_200 = [0]

#parasolic motion detectors complete
#ps_all_file = open('data/mo_det_cal/p_max_spikes_all.txt','a+')

#store_sp=[]
#store_sp_200=[]

#firstq = True
#for q in range(1,31):
#    s = tp.FindNearestElement(out_p_r_60,[7.5,4*float(q/2.)])
#    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
#    s_gids += [s]
    
#    dSD = nest.GetStatus(s,keys="events")[0]
#    evs = dSD["senders"]
#    ts = dSD["times"]
#    #pyl.figure(2)
#    #pyl.plot(ts, evs, ".")
#    #pyl.show()
    
#    if evs.any():
#        #ev = nest.GetStatus(mult)[0]['events']
#        #t = ev['times']
#        #r = ev['rate']

#        sp = nest.GetStatus(s)[0]['events']['times']
#        #plt.subplot(221)
#        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
#        store_sp +=[h]
#        h200=0
#        for spike in sp:
#            if spike>200.:
#                h200+=1
#        store_sp_200 +=[h200]
#        qa = [4*q for i in range(h)]
#        if firstq == True:
#            ax.plot(qa,sp,'k*',label='general motion detectors')
#            firstq = False
#        else:
#            ax.plot(qa,sp,'k*',label='_nolegend_')
#        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
#        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
#        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
#        #plt.plot(t, r, color='b')
#        #plt.step(e[:-1], h , color='b', where='post')
#        #plt.title('PST histogram and firing rates')
#        #plt.ylabel('Spikes per second')

#        #plt.subplot(223)
#        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
#                    #histtype='step', color='b')
#        #plt.title('ISI histogram')
#        #plt.show()
#    else:
#        store_sp +=[0]
#        store_sp_200 +=[0]

#ps_all_file.close()
##plt.show()
#max_spikes_ud = max(store_sp)
#if max_spikes_ud == 0:
#    max_spikes_ud = [0]
#print max_spikes_ud

#max_spikes_ud_200 = [max(store_sp_200)]
#if max_spikes_ud_200 == 0:
#    max_spikes_ud_200 = [0]
##print max_spikes_200

#up-down-part-120--------------------------------------------------------------------------------------------------------------------

#parasolic motion detectors up 120
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_up_120_'+str(handle_name)+'.txt','w+')

firstq = True
#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_120_up,[pos[0],pos[1]])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='120 downward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_up_120 = max(store_sp)
if max_spikes_up_120 == 0:
    max_spikes_up_120 = [0]
print max_spikes_up_120

max_spikes_up_120_200 = [max(store_sp_200)]
if max_spikes_up_120_200 == 0:
    max_spikes_up_120_200 = [0]

#parasolic motion detectors down 120
store_sp=[]
store_sp_200=[]

directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
if not os.path.exists(directory):
    os.makedirs(directory)
this_file = open(str(directory)+'/parasols_md_down_120_'+str(handle_name)+'.txt','w+')

firstq = True

#for q in range(1,31):
for pos in gp_pos:
    s = tp.FindNearestElement(out_p_r_120_down,[pos[0],pos[1]]) 
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        this_file.write(str(sp)+'\n')
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        all_movement_spikes+=h
        h200=0
        for spike in sp:
            if spike>200.:
                h200+=1
        store_sp_200 +=[h200]
        all_movement_spikes_200+=h200
        
        #qa = [4*q for i in range(h)]
        #if firstq == True:
        #    ax.plot(qa,sp,'rx',label='120 downward-motion detectors')
        #    firstq = False
        #else:
        #    ax.plot(qa,sp,'rx',label='_nolegend_')
    else:
        store_sp +=[0]
        store_sp_200 +=[0]
        this_file.write('[0]\n')

this_file.close() 

max_spikes_down_120 = max(store_sp)
if max_spikes_down_120 == 0:
    max_spikes_down_120 = [0]
print max_spikes_down_120

max_spikes_down_120_200 = [max(store_sp_200)]
if max_spikes_down_120_200 == 0:
    max_spikes_down_120_200 = [0]



#image/output-part---------------------------------------------------------------------------------------------------------------

#print weight,delay

#plt.title('midget spikes for motion detectors')
#plt.xlabel('position of neuron in arcmin')
#plt.ylabel('time $t$ in ms')
#plt.legend()
'''
plt.xlim([0,120])
plt.ylim([0,t_end-t_start])

#ax.set_title('spiking times uniformly distributed with $\sigma_{I_e} = 10$pA')
ax.set_xlabel('position of neuron in arcmin')
#ax.get_xaxis().set_visible(False)
ax.set_ylabel('time $t$ in ms')
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

#ax.legend(loc=4, fancybox=False, fontsize=8, shadow=False)
# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), 
          fancybox=False, fontsize=8, shadow=False, ncol=2)

plt.savefig('img/mo_det_cal/x_t_parasols'+handle_name+'_90deg.pdf')
plt.savefig('img/mo_det_cal/x_t_parasols_'+handle_name+'_90deg.pgf')
out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '_x_t_parasols'+ str(handle_name) + '.pdf'
plt.savefig(out)
plt.show()
plt.close()
'''

#and maximum spikes for one cell
max_spikes_total = max(max_spikes_left[0],max_spikes_right[0],max_spikes_up_60[0],max_spikes_down_60[0],max_spikes_up_120[0],max_spikes_down_120[0])
max_spikes_total_200 = max(max_spikes_left_200[0],max_spikes_right_200[0],max_spikes_up_60_200[0],max_spikes_down_60_200[0],max_spikes_up_120_200[0],max_spikes_down_120_200[0])
print max_spikes_total

ps_file = open('data/poletti2010/m_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(max_spikes[0])+'\t'+str(m_in_max_spikes[0])+'\t'+ str(max_spikes_left[0])+'\t'+ str(max_spikes_right[0])+'\t'+ str(max_spikes_200[0])+'\t'+ str(m_in_max_spikes_200[0])+'\t'+ str(max_spikes_left_200[0])+'\t'+ str(max_spikes_right_200[0])+'\n')
#ps_file = open('data/mo_det_cal/m_max_spikes.txt','a+')
ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(p_in_max_spikes[0])+'\t'+ str(max_spikes_left[0])+'\t'+ str(max_spikes_right[0])+'\t'+ str(max_spikes_up_60[0])+'\t'+ str(max_spikes_down_60[0])+'\t'+ str(max_spikes_up_120[0])+'\t'+ str(max_spikes_down_120[0])+'\t'+str(max_spikes_total)+'\t'+ str(all_movement_spikes)+'\n')
ps_file.close()

ps_file_200 = open('data/poletti2010/m_max_spikes_200.txt','a+')
#ps_file_200 = open('data/mo_det_cal/m_max_spikes_200.txt','a+')
ps_file_200.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+ str(p_in_max_spikes_200[0])+'\t'+ str(max_spikes_left_200[0])+'\t'+ str(max_spikes_right_200[0])+'\t'+ str(max_spikes_up_60_200[0])+'\t'+ str(max_spikes_down_60_200[0])+'\t'+ str(max_spikes_up_120_200[0])+'\t'+ str(max_spikes_down_120_200[0])+'\t'+str(max_spikes_total_200)+'\t'+ str(all_movement_spikes_200)+'\n')
#ps_file_200.write(str(weight)+'\t'+str(delay)+'\t'+ str(vel)+'\t'+ str(p_in_max_spikes_200[0])+'\t'+ str(max_spikes_left_200[0])+'\t'+ str(max_spikes_right_200[0])+'\t'+ str(max_spikes_up_60_200[0])+'\t'+ str(max_spikes_down_60_200[0])+'\t'+ str(max_spikes_up_120_200[0])+'\t'+ str(max_spikes_down_120_200[0])+'\t'+str(max_spikes_total_200)+'\t'+ str(all_movement_spikes_200)+'\n')
ps_file_200.close()
