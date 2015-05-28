import sys
import math

if (len(sys.argv) != 7):
  print("usage:  gen.py #x_cells #y_cells #z_cells base#particlespercell beamradius cutoff ")


# input number of cells and approx number of particles
nx = int(sys.argv[1])  # num cells in x
ny = int(sys.argv[2])  # num cells in y
nz = int(sys.argv[3])  # num cells in z
np = int(sys.argv[4])  # base number of particles per cell to use
cutradius = float(sys.argv[5])  # radius of the beam; no particles outside
cutz = float(sys.argv[6])       # z value after which there are no particles

totalcells = nx * ny * nz

# fixed box dimensions
# math below assumes x, y and z are centered at zero
minx = -2.
maxx =  2.
miny = -2.
maxy =  2.
minz =  0.
maxz =  10.

# Center of circle defining beam input
cx = 0.5 * (maxx + minx)
cy = 0.5 * (maxy + miny)

# deltas
dx = (maxx - minx)
dy = (maxy - miny)
dz = (maxz - minz)

# mesh spacing
hx = dx / nx
hy = dy / ny
hz = dz / nz

halfhx = hx * 0.5
halfhy = hy * 0.5
halfhz = hz * 0.5

# open files for coordinates and graph
fc = open(str(nx)+"x"+str(ny)+"x"+str(nz)+".coords", "w")
fg = open(str(nx)+"x"+str(ny)+"x"+str(nz)+".graph", "w")
fs = open(str(nx)+"x"+str(ny)+"x"+str(nz)+".stats", "w")
WEIGHTCODE = " 0 010 2"  # 0 edges; 010 for vtx wghts; 2 for two wghts per vtx
fg.write(str(totalcells) + WEIGHTCODE + '\n')

fs.write(sys.argv[0]+" "+ sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]+" "+sys.argv[4]+" "+sys.argv[5]+" "+sys.argv[6]+"\n\n")
fs.write("nx,ny,nz = (" +str(nx)+","+str(ny)+","+str(nz)+")\n")
fs.write("dx,dy,dz = (" +str(dx)+","+str(dy)+","+str(dz)+")\n")
fs.write("minx,miny,minz = (" +str(minx)+","+str(miny)+","+str(minz)+")\n")
fs.write("maxx,maxy,maxz = (" +str(maxx)+","+str(maxy)+","+str(maxz)+")\n")
fs.write("base np="+str(np)+"\n")
fs.write("cutradius="+str(cutradius)+"\n")
fs.write("cutz="+str(cutz)+"\n\n\n")

# generate mesh
totalparticles = 0
maxwgt = 0
minwgt = 1234567
maxcnt = 0
mincnt = 0
noutsideradius = 0
npastzcutoff = 0

for ix in range(nx):
    x = minx + ix * hx + halfhx

    for iy in range(ny):
        y = miny + iy * hy + halfhy

        for iz in range(nz):
            z = minz + iz * hz + halfhz
            if (math.sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) > cutradius):
                # cell is outside beam radius in x-y plane
                wgt = 0 
                noutsideradius += 1
            else:
                # cell is inside beam radius in x-y plane
                distz = z - minz
                if (distz > (minz + cutz)):
                    # cell is beyond the z-distance cutoff
                    wgt = 0
                    npastzcutoff += 1
                else:
                    if (distz > 0):
                        wgt = int(np / distz)
                    else:
                        wgt = np

            totalparticles += wgt

            if (wgt > maxwgt):
                maxwgt = wgt
                maxcnt = 0
            if (wgt == maxwgt):
                maxcnt += 1
            if (wgt < minwgt):
                minwgt = wgt
                mincnt = 0
            if (wgt == minwgt):
                mincnt += 1
            fc.write("%1.4f "%x + "%1.4f "%y + "%1.4f\n"%z)
            fg.write(str(wgt) + ' 1\n')

fs.write("Total: " + str(totalcells) + " cells\n")
fs.write("       " + str(totalparticles) + " particles\n")
fs.write("Max part/cell " + str(maxwgt) + " in " + str(maxcnt) + " cells\n")
fs.write("Min part/cell " + str(minwgt) + " in " + str(mincnt) + " cells\n")
fs.write("Avg part/cell " + str(totalparticles / totalcells) + "\n")
fs.write("Num cells outside beam radius " + str(noutsideradius) + "\n")
fs.write("Num cells past z cutoff       " + str(npastzcutoff) + "\n")

fc.close()
fg.close()
fs.close()
