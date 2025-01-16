#########Function to define the pair of bonds between two side-by-side lattices
def bondPairs(Nx,Ny):
  N = Nx*Ny
  bonds = []

  for i in range(1,N):
    if i%Ny == 1 and i<N-Ny:
        p1 = [i,i+1]
        bonds.append(p1)
        p1 = [i,i+Ny]
        bonds.append(p1)
        p1 = [i,i+Ny+1]
        bonds.append(p1)
        p1 = [i,i+Ny-1]
        bonds.append(p1)
        p1 = [i,i+2*Ny-1]
        bonds.append(p1)

    elif i%Ny == 0 and i <N:
        p1 = [i,i+Ny]
        bonds.append(p1)

    elif i>N-Ny:
        p1 = [i,i+1]
        bonds.append(p1)
        if i==N-Ny+1:
            p1 = [i,i+Ny-1]
            bonds.append(p1)

    elif i%2 == 1:
        p1 = [i,i+1]
        bonds.append(p1)
        p1 = [i,i+Ny-1]
        bonds.append(p1)
        p1 = [i,i+Ny]
        bonds.append(p1)
        p1 = [i,i+Ny+1]
        bonds.append(p1)

    elif i%2 == 0:
        p1 = [i,i+1]
        bonds.append(p1)
        p1 = [i,i+Ny]
        bonds.append(p1)

  return bonds

############################# Function to create reciprocal coordinate
def reciprocal_lattice(a1,a2,a3):
  reciproc_coord = []
  V = np.dot(a1,np.cross(a2,a3))
  b1 = np.cross(a2,a3)
  b2 = np.cross(a3,a1)
  b3 = np.cross(a1,a2)

  for b in[b1,b2,b3]:
    reciproc_coord.append(b*2*np.pi/V)
  return reciproc_coord

############################# Generating the triangular coordinate
def triangular_coord(Nx,Ny,a1,a2):
  sites = 1 #auxiliary variables to denote the sites
  x = 1 #initial position
  y = 0 #initial position
  site_coord = {}
  for i in range(0,Nx):
    x_art = x
    y = 0
    for j in range(0,Ny):
      if j != 0:
        if (j+1)%2 != 0:
          x_art += a2[0]
        elif ((j+1)%2 == 0):
          x_art -= a2[0]
      site_coord[sites] = [x_art,y]
      y += np.round(a2[1],3)
      sites += 1
    x += a1[0]
  return site_coord

############################# Function to calculate the correlator
def correlator_momentum(kx,ky,site_coord,C):
  C_k = 0
  for i in range(0,N):
    for j in range(0,N):
      x_ij = (site_coord[bnd[0]][0]-site_coord[bnd[1]][0])
      y_ij = (site_coord[bnd[0]][1]-site_coord[bnd[1]][1])
      if site_coord[i+1] == site_coord[j+1]:
        C_k += 1/N * C[i][j]
      else:
        C_k += 1/N * np.exp(kx*x_ij + ky*y_ij) * C[i][j]
  return C_k
