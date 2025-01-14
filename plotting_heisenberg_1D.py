import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import matplotlib as mpl

##Generating lattice for testing
N = 10

fig,ax = plt.subplots()

for i in range(0,N):
  ax.scatter(i+1,5,50,color="k")
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k")

##Generating lattice for testing Hamiltonian results
N = 40

df_corr = pd.read_csv("...csv") ##use the dataset that have  SzSz, S+S- and S-S+
df_corr["<SiSj>"] = df_corr.loc[:,"<SzSz>"] + 1/2 * (df_corr.loc[:,"<S+S->"] + df_corr.loc[:,"<S-S+>"]) ##Calcualte SiSj

C_szsz = np.zeros([N,N])
C_spsm = np.zeros([N,N])
C_smsp = np.zeros([N,N])
C_sisj = np.zeros([N,N])

k = 0

##Placing the results to the coordinate of the maps
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

################### Plotting the results of S
fig,ax = plt.subplots()

for i in range(0,N):
  ax.scatter(i+1,5, df_sz["<Sz>"][i]*1E7, color="b", zorder=2)
  ax.scatter(i+1,5, -df_sz["<Sz>"][i]*1E7, color="r", zorder=2)
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k", linewidth=1, zorder=1)


ax.axes.get_yaxis().set_visible(False)
plt.xlabel("x (in a)")
plt.xticks( range(0,10,1) )

plt.title("Magnetization")
plt.savefig('...png', dpi=1500) ##Put name of file
plt.show()

######################################## Plot additional graphs for SiSj separately
ind_foc = 5
fig,ax = plt.subplots()

for i in range(0,N):
  if i != ind_foc-1:
    ax.scatter(i+1,5, C_szsz[ind_foc-1][i]*500, color="b", zorder=2)
    ax.scatter(i+1,5, -C_szsz[ind_foc-1][i]*500, color="r", zorder=2)
  else:
    ax.scatter(i+1,5, 100, color = "g",zorder=2)
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k", linewidth=1, zorder=1)


ax.axes.get_yaxis().set_visible(False)
plt.xlabel("x (in a)")
plt.xticks( range(0,10,1) )

plt.title("<SzSz>")
#plt.savefig('/content/SzSz_heisenberg.png', dpi=2000)
plt.show()

fig,ax = plt.subplots()

for i in range(0,N):
  if i != ind_foc-1:
    ax.scatter(i+1,5, C_spsm[ind_foc-1][i]*500, color="b", zorder=2)
    ax.scatter(i+1,5, -C_spsm[ind_foc-1][i]*500, color="r", zorder=2)
  else:
    ax.scatter(i+1,5, 100, color = "g",zorder=2)
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k", linewidth=1, zorder=1)


ax.axes.get_yaxis().set_visible(False)
plt.xlabel("x (in a)")
plt.xticks( range(0,11,1) )

plt.title("<S+S->")
#plt.savefig('/content/SpSm_heisenberg.png', dpi=2000)
plt.show()

fig,ax = plt.subplots()

for i in range(0,N):
  if i != ind_foc-1:
    ax.scatter(i+1,5, C_smsp[ind_foc-1][i]*500, color="b", zorder=2)
    ax.scatter(i+1,5, -C_smsp[ind_foc-1][i]*500, color="r", zorder=2)
  else:
    ax.scatter(i+1,5, 100, color = "g",zorder=2)
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k", linewidth=1, zorder=1)


ax.axes.get_yaxis().set_visible(False)
plt.xlabel("x (in a)")
plt.xticks( range(0,11,1) )

plt.title("<S-S+>")
#plt.savefig('/content/SmSp_heisenberg.png', dpi=2000)
plt.show()

fig,ax = plt.subplots()

for i in range(0,N):
  if i != ind_foc-1:
    ax.scatter(i+1,5, C_sisj[ind_foc-1][i]*500, color="b", zorder=2)
    ax.scatter(i+1,5, -C_sisj[ind_foc-1][i]*500, color="r", zorder=2)
  else:
    ax.scatter(i+1,5, 100, color = "g",zorder=2)
  if i != 0:
    ax.plot([i,i+1],[5,5], color = "k", linewidth=1, zorder=1)


ax.axes.get_yaxis().set_visible(False)
plt.xlabel("x (in a)")
plt.xticks( range(0,11,1) )

plt.title("<SiSj>")
#plt.savefig('/content/SiSj_heisenberg.png', dpi=1500)
plt.show()


##########################Plotting energy per site testing
energy_per_site = [-0.375,-0.4040063509461097,-0.4155961889813205,-0.4218665748357683,-0.42580351833594304,-0.4285075496582354,-0.43048032530248886,-0.43198356295319895,-0.43316726883297146,-0.43412365256161456]
x = [2,4,6,8,10,12,14,16,18,20]

plt.figure()
plt.xticks( range(0,21,2) )
plt.ylabel("Energy per site")
plt.xlabel("Number of site")
plt.title("Energy per site")

plt.scatter(x,energy_per_site)
plt.savefig('/content/heisenberg_energy_per_site.png', dpi=1500)


plt.figure()
plt.xticks( range(0,21,2) )
plt.ylim([-1,0])
plt.ylabel("Energy per site")
plt.xlabel("Number of site")
plt.title("Energy per site")

plt.scatter(x,energy_per_site)
plt.savefig('/content/heisenberg_energy_per_site_larger.png', dpi=1500)
