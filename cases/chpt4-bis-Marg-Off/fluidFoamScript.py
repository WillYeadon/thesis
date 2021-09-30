from fluidfoam import readscalar
from fluidfoam.readpostpro import readforce



sol = '.'
#timename = '0'

force = readforce(sol, time_name = 'mergeTime')
#output = readscalar(sol, timename, 'alpha.Phase1')

import matplotlib.pyplot as plt

plt.figure()

plt.plot(force[:, 0], force[:, 1])

# Setting axis labels
plt.xlabel('t (s)')
plt.ylabel('p (Pa)')

# add grid
plt.grid()


#import matplotlib.pyplot as plt
#plt.figure()
#plt.plot(output)
#plt.savefig('foo.png')