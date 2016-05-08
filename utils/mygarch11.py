#############################################################################
#
# Example of a simple garch 1 1 model
#
# Copyright (C) 2016 Yingfeng Yu < yuyingfeng (at) cueb.edu.cn >
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
##############################################################################

# Import some libraries
import Quandl
import matplotlib.pylab          as plt
import numpy                     as np

# Set the random seed to replicate results in tutorial
np.random.seed( 10 );



##############################################################################
# Load data
##############################################################################
d = Quandl.get("NASDAQOMX/OMXS30", trim_start="2012-01-02", trim_end="2014-01-02")
#d = qu.get("NASDAQOMX/OMXS30", trim_start="2012-01-02", trim_end="2014-01-02")
print("============================PASS Quandl.get=================================================")
y = 100 * np.diff(np.log(d['Index Value']))


##############################################################################
# Define the model
##############################################################################
nIter=len(y);
print(nIter)

alpha=0.87;
beta=0.1;
vol   = np.zeros((nIter,1));

#mu=np.mean(y)
mu=1-alpha-beta;
vol[0]=np.std(y)+mu;

for tt in range(1, nIter):

        tmp=mu+alpha*(vol[tt-1]**2)+beta*(y[tt-1]**2);
        vol[tt]=np.sqrt(tmp);


plt.figure(1);

# Plot the log-returns
plt.subplot(2,1,1);
plt.plot(y,color = '#1B9E77', linewidth=1.5);
plt.xlabel("time"); plt.ylabel("log-return")

# Plot the volatility
plt.subplot(2,1,2);
plt.plot(vol,color = '#E7298A', linewidth=1.5);
plt.xlabel("time"); plt.ylabel("volatility")
plt.show()
