import matplotlib.pyplot as plt
import csv
import numpy as np

title = "stats for PRICE, for changing M"
file_name = "price.csv"

est = [];
time = [];
err_mean = [];
err_var = [];
M = [];

f = open(file_name, "r")
reader = csv.reader(f)
for row in reader:
    est.append( float( row[0] ) )
    time.append( float( row[1] ) )
    err_mean.append( float( row[2] ) )
    err_var.append( float( row[3] ) )
    M.append( float( row[4] ) )

plt.figure()
plt.title(title)
plt.plot( M,est )
plt.xlabel("number of simulations")
plt.ylabel("Estimate")
plt.show()

plt.figure()
plt.title(title)
plt.plot( M,time )
plt.xlabel("number of simulations")
plt.ylabel("Time")
plt.show()

plt.figure()
plt.title(title)
plt.plot( M,err_mean )
plt.xlabel("number of simulations")
plt.ylabel("Err mean")
plt.show()

plt.figure()
plt.title(title)
plt.plot( M,err_var )
plt.xlabel("number of simulations")
plt.ylabel("Err var")
plt.show()