import itertools
import numpy as np



hello = [[0.5,0.8,0.6],[0.3,0.5,0.1]]
u = [3,4,2,4,5,8,2,9,0]
'''
for i,p in zip(enumerate(hello)):
    #r= p*q
    print p,i
'''  
depos=[(0.,0.866),(.5,0.),(1.5,0.),(2.,.866)]
for i in zip(u,itertools.cycle(depos)):
    print i[1][0]
    
'''
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

def generate_data():
    a = np.arange(25).reshape(5, 5)
    b = 10 * np.random.rand(5, 5)
    return a - b 

def update(data):
    mat.set_data(data)
    return mat 

def data_gen():
    while True:
        yield generate_data()

fig, ax = plt.subplots()
mat = ax.matshow(generate_data())
plt.colorbar(mat)
ani = animation.FuncAnimation(fig, update, data_gen, interval=500,
                              save_count=50)
plt.show()

ani.save('animation.mp4')

ani.save('animation.mp4', clear_temp=False)

convert *.png animation.gif