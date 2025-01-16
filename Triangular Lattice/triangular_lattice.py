import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import matplotlib as mpl

######## Testing Coordinate
a_const = 1 #lattice constant
a = [[a_const,0,0],[a_const/2,(a_const*np.sqrt(3)/2),0],[0,0,1]]
b = reciprocal_lattice(a[0],a[1],a[2])
Nx = 10
Ny = 4
N = Nx * Ny
site_coord = triangular_coord(Nx,Ny,a[0],a[1])
bond = bondPairs(Nx,Ny)

fig, ax = plt.subplots()

for i in range(1,Nx*Ny+1):
  ax.scatter(site_coord[i][0],site_coord[i][1], 5, color="k")
  ax.annotate(i, (site_coord[i][0],site_coord[i][1]))
for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const : 
    ax.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k")

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")
#plt.savefig('/content/triangular_geometry.png', dpi=1200)
plt.show()

#################################
df_corr = pd.read_csv("/content/test5.csv")
df_corr["<SiSj>"] = df_corr.loc[:,"<SzSz>"] + 1/2 * (df_corr.loc[:,"<S+S->"] + df_corr.loc[:,"<S-S+>"])

C_szsz = np.zeros([N,N])
C_spsm = np.zeros([N,N])
C_smsp = np.zeros([N,N])
C_sisj = np.zeros([N,N])

k = 0

for i in range(0,N):
  for j in range(i,N):
    C_szsz[i][j] = df_corr["<SzSz>"].iloc[k]
    C_szsz[j][i] = C_szsz[i][j]

    C_spsm[i][j] = df_corr["<S+S->"].iloc[k]
    C_spsm[j][i] = C_spsm[i][j]

    C_smsp[i][j] = df_corr["<S-S+>"].iloc[k]
    C_smsp[j][i] = C_smsp[i][j]

    C_sisj[i][j] = df_corr["<SiSj>"].iloc[k]
    C_sisj[j][i] = C_sisj[i][j]

    k += 1

S = 0
for i in range(0,N):
  for j in range(0,N):
    S += C_sisj[i][j]

print(S)

#df_sz = pd.read_csv("/content/magnetization.csv")


##Printing the graph coordinate
plt.figure()

for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const :
    plt.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k", linewidth=1, zorder=1)


for i in range(1,25):
  plt.scatter(site_coord[i][0],site_coord[i][1], df_sz["<Sz>"][i-1]*500, color="b", zorder=2)
  plt.scatter(site_coord[i][0],site_coord[i][1], -df_sz["<Sz>"][i-1]*500, color="r", zorder=3)

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")

plt.title("Magnetization")
#plt.savefig('/content/Magnetization.png', dpi=1200)
plt.show()


#### Printing the Results on the Triangular Lattice - SzSz
ind_foc = 18
plt.figure()

for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const :
    plt.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k", linewidth=1, zorder=1)

for i in range(1,41):
  if i != ind_foc:
    plt.scatter(site_coord[i][0],site_coord[i][1], C_szsz[ind_foc-1][i-1]*500, color="b", zorder=2)
    plt.scatter(site_coord[i][0],site_coord[i][1], -C_szsz[ind_foc-1][i-1]*500, color="r", zorder=3)
  else:
    plt.scatter(site_coord[i][0],site_coord[i][1], 100, color = "g",zorder=2)

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")

plt.title("<SzSz>")
plt.show()

plt.figure()
for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const :
    plt.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k", linewidth=1, zorder=1)

for i in range(1,41):
  if i != ind_foc:
    plt.scatter(site_coord[i][0],site_coord[i][1], C_spsm[ind_foc-1][i-1]*500, color="b", zorder=2)
    plt.scatter(site_coord[i][0],site_coord[i][1], -C_spsm[ind_foc-1][i-1]*500, color="r", zorder=3)
  else:
    plt.scatter(site_coord[i][0],site_coord[i][1], 100, color = "g",zorder=2)

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")

plt.title("<S+S->")
plt.show()

plt.figure()
for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const :
    plt.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k", linewidth=1, zorder=1)

for i in range(1,41):
  if i != ind_foc:
    plt.scatter(site_coord[i][0],site_coord[i][1], C_smsp[ind_foc-1][i-1]*5000, color="b", zorder=2)
    plt.scatter(site_coord[i][0],site_coord[i][1], -C_smsp[ind_foc-1][i-1]*5000, color="r", zorder=3)
  else:
    plt.scatter(site_coord[i][0],site_coord[i][1], 100, color = "g",zorder=2)

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")

plt.title("<S-S+>")
plt.show()

plt.figure()
for bnd in bond:
  if np.round(np.sqrt((site_coord[bnd[0]][0]-site_coord[bnd[1]][0])**2 + (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])**2)) == 1*a_const :
    plt.plot([site_coord[bnd[0]][0],site_coord[bnd[1]][0]],[site_coord[bnd[0]][1],site_coord[bnd[1]][1]], color="k", linewidth=1, zorder=1)


for i in range(1,41):
  if i != ind_foc:
    plt.scatter(site_coord[i][0],site_coord[i][1], C_sisj[ind_foc-1][i-1]*5000, color="b", zorder=2)
    plt.scatter(site_coord[i][0],site_coord[i][1], -C_sisj[ind_foc-1][i-1]*5000, color="r", zorder=3)
  else:
    plt.scatter(site_coord[i][0],site_coord[i][1], 100, color = "g",zorder=2)

plt.xlabel("x (in a)")
plt.ylabel("y (in a)")
#plt.savefig('/content/SiSj_11.png', dpi=1200)
plt.title("<SiSj>")
plt.show()
