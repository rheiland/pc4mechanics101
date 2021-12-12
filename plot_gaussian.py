import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# define constants
mu = 998.8 
sigma = 73.10
x1 = 900
x2 = 1100

# calculate Z-transform of the lower and upper bound using the mean and standard deviation defined above.
# calculate the z-transform
z1 = ( x1 - mu ) / sigma
z2 = ( x2 - mu ) / sigma

x = np.arange(z1, z2, 0.001) # range of x in spec
#x_all = np.arange(-10, 10, 0.001) # entire range of x, both in and out of spec
xmin = -30
xmax = -xmin
# x_all = np.arange(-10, 10, 0.001) # entire range of x, both in and out of spec
x_all = np.arange(xmin, xmax, 0.1) # entire range of x, both in and out of spec
# mean = 0, stddev = 1, since Z-transform was calculated
y = norm.pdf(x,0,1)
y2 = norm.pdf(x_all,0,1)
y2 *= 5

# build the plot
fig, ax = plt.subplots(figsize=(8,4))
#plt.style.use('fivethirtyeight')
ax.plot(x_all,y2)

# ax.fill_between(x,y,0, alpha=0.3, color='b')
# ax.fill_between(x_all,y2,0, alpha=0.1)
# ax.set_xlim([-4,4])
ax.set_xlim([xmin,xmax])
# ax.set_xlabel('# of Standard Deviations Outside the Mean')
# ax.set_yticklabels([])
ax.set_title('Normal Gaussian Curve')
plt.show()
