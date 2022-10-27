import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import mpl_toolkits.mplot3d.art3d as art3d
import uproot

# Open data file
fIn = uproot.open('S2-tuple.root')

# Select trees of interest
clusterTree = fIn['DefaultReconstruction/clusters']
muonTree = fIn['DefaultReconstruction/muons']

# Select leaves of interest
clusterTreeEntries=clusterTree.arrays(['start', 'rejected', 'MuonTag', 'NSTag', 'ESTag', 'MuonTreePos', 'chanID', 'chanActive', 'cubes', 'wfTreePos', 'chanAmplitude', 'onionLayer'])
muonTreeEntries = muonTree.arrays(['htime','startCube','endCube','type','stopped', 'mx', 'cx', 'my', 'cy'])


# --------------------------------------------------------------------------------------------
# set scales of axes:

# Make sure these are floating point values:
scale_x = 5.0
scale_y = 1.6
scale_z = 1.6

# Axes are scaled down to fit in scene
max_scale = max(scale_x, scale_y, scale_z)

scale_x = scale_x/max_scale
scale_y = scale_y/max_scale
scale_z = scale_z/max_scale

# Create scaling matrix
scale = np.array([[scale_x, 0, 0, 0], [0, scale_y, 0, 0], [0, 0, scale_z, 0], [0, 0, 0, 0.6]])

def get_proj_scale(self):
    """
    Create the projection matrix from the current viewing position.
    elev stores the elevation angle in the z plane
    azim stores the azimuth angle in the x,y plane
    dist is the distance of the eye viewing point from the object point.
    """
    relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

    xmin, xmax = self.get_xlim3d()
    ymin, ymax = self.get_ylim3d()
    zmin, zmax = self.get_zlim3d()

    # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
    worldM = proj3d.world_transformation(xmin, xmax, ymin, ymax, zmin, zmax)

    # look into the middle of the new coordinates
    R = np.array([1, 0.4, -0.3])

    xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
    yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
    zp = R[2] + np.sin(relev) * self.dist
    E = np.array((xp, yp, zp))

    self.eye = E
    self.vvec = R - E
    self.vvec = self.vvec / proj3d.mod(self.vvec)

    if abs(relev) > np.pi/2:
        # upside down
        V = np.array((0, 0, -1))
    else:
        V = np.array((0, 0, 1))
    zfront, zback = -self.dist, self.dist

    viewM = proj3d.view_transformation(E, R, V)
    perspM = proj3d.persp_transformation(zfront, zback)
    M0 = np.dot(viewM, worldM)
    M = np.dot(perspM, M0)

    return np.dot(M, scale)
    
# -------------------------------------------------------------------------------------------
# determine channel
def chanIDTranslation(iC):
    z = int(iC/64)
    iC = iC % 64

    if (iC < 8):
        horNotVer = False
        x = 7-iC
        y = 16

    elif (iC < 16):
        horNotVer = True
        x = -1
        y = 15-(iC-8)
        if (iC == 14): 
            y = 8
        if (iC == 15):
            y = 9

    elif (iC < 24):
        horNotVer = True
        x = 16
        y = iC-8

    elif (iC < 32):
        horNotVer = False
        x = 8 + (31-iC)
        y = 16

    elif (iC < 40):
        horNotVer = False
        x = iC-24
        y = -1

    elif (iC < 48):
        horNotVer = True
        x = 16
        y = 8-(48-iC)
        if (iC == 47):
            y=6
        if (iC == 46): 
            y=7

    elif (iC < 56):
        horNotVer = True
        x = -1
        y = 8+(47-iC)

    else:
        horNotVer = False
        x = iC-56
        y = -1

    return x, y, z, horNotVer
    
# --------------------------------------------------------------------------------------------
def VisualFromFile(filename, eventNumber, figName, eventOneName, eventTwoName=None, isMuon=False, numberOfEventsTogether=2, alpha=0.8):
    file = open(filename, 'r')

    # read in chosen event from file
    i = 0
    for line in file:
        if i == eventNumber:
            chosen = line
        i += 1
    file.close()
    if numberOfEventsTogether == 1:
        iEntry1 = int(chosen)
    elif numberOfEventsTogether == 2:
        events = chosen.split('\t')
        iEntry1 = int(events[0])
        iEntry2 = int(events[1])
    # --------------------------------------------------
    # plotten
    Axes3D.get_proj = get_proj_scale

    fig = plt.figure(figsize=(9, 7))
    gs = gridspec.GridSpec(nrows=2, ncols=1)

    ax = fig.add_subplot(gs[0:, 0], projection='3d')
    ax1 = fig.add_subplot(gs[1, 0])

    # normalisation for color map
    nrm = mpl.colors.Normalize(0, 16000)
    cmap = 'viridis'
    colorsMuon = cm.viridis(nrm(clusterTreeEntries['chanAmplitude'][iEntry1]))
    colorsElectron = cm.viridis(nrm(clusterTreeEntries['chanAmplitude'][iEntry2]))

    matrix1 = [[ [0 for col in range(16)] for col in range(16)] for row in range(50)]
    z = np.zeros(len(clusterTreeEntries['cubes'][iEntry1]))

    for i in range(len(clusterTreeEntries['cubes'][iEntry1])):
        matrix1[int(clusterTreeEntries['cubes'][iEntry1][i][2])][int(clusterTreeEntries['cubes'][iEntry1][i][0])][int(clusterTreeEntries['cubes'][iEntry1][i][1])] = 1
        z[i] = clusterTreeEntries['cubes'][iEntry1][i][2]
    if isMuon:
        Zspace = np.linspace(min(z)-5, max(z)+5, 1000)
        jEntry1 = clusterTreeEntries['MuonTreePos'][iEntry1]

        X = muonTreeEntries['mx'][jEntry1] * Zspace + muonTreeEntries['cx'][jEntry1]+0.5
        Y = muonTreeEntries['my'][jEntry1] * Zspace + muonTreeEntries['cy'][jEntry1]+0.5

        ax.plot3D(Zspace, X, Y, linewidth=0.7, color='#0000FF', label='{}'.format(eventOneName))

    """
    Note: using HTML colors allows us to set the transparency by adding the last two digits (the alpha value)
    """
    ax.voxels(np.array(matrix1), facecolors='#0000FF75', edgecolor='#08298A75')

    for n, iC in enumerate(clusterTreeEntries['chanID'][iEntry1]):
        x, y, z, horNotVer = chanIDTranslation(iC)
        if horNotVer and (x == 16):
            Xc = np.linspace(-1, 16, 10)
            Yc = np.repeat(y + 1, 10)
            Zc = np.repeat(z + 1, 10)

            rectangle = plt.Rectangle((z, y), 1, 1, fc=colorsMuon[n], ec=colorsMuon[n], alpha=alpha)  # projection
            ax.add_patch(rectangle)
            art3d.pathpatch_2d_to_3d(rectangle, z=16, zdir="y")

        elif horNotVer and (x == -1):
            Xc = np.linspace(-1, 16, 10)
            Yc = np.repeat(y, 10)
            Zc = np.repeat(z, 10)

            rectangle = plt.Rectangle((z, y), 1, 1, fc=colorsMuon[n], ec=colorsMuon[n], alpha=alpha)  # projection
            ax.add_patch(rectangle)
            art3d.pathpatch_2d_to_3d(rectangle, z=-1, zdir="y")

        elif (not horNotVer) and (y == 16):
            Xc = np.repeat(x + 1, 10)
            Yc = np.linspace(-1, 16, 10)
            Zc = np.repeat(z + 1, 10)

            rectangle = plt.Rectangle((z, x), 1, 1, fc=colorsMuon[n], ec=colorsMuon[n], alpha=alpha)  # projection
            ax.add_patch(rectangle)
            art3d.pathpatch_2d_to_3d(rectangle, z=16, zdir="z")
        else:
            Xc = np.repeat(x, 10)
            Yc = np.linspace(-1, 16, 10)
            Zc = np.repeat(z, 10)

            rectangle = plt.Rectangle((z, x), 1, 1, fc=colorsMuon[n], ec=colorsMuon[n], alpha=alpha)  # projection
            ax.add_patch(rectangle)
            art3d.pathpatch_2d_to_3d(rectangle, z=-1, zdir="z")

        ax.plot3D(Zc, Xc, Yc, linewidth=0.13, color='#0000FF')  # channels

    # if there are 2 events to be plotted together
    if numberOfEventsTogether == 2:
        matrix2 = [[ [0 for col in range(16)] for col in range(16)] for row in range(50)]
        z = np.zeros(len(clusterTreeEntries['cubes'][iEntry2]))

        for j in range(len(clusterTreeEntries['cubes'][iEntry2])):
            matrix2[int(clusterTreeEntries['cubes'][iEntry2][j][2])][int(clusterTreeEntries['cubes'][iEntry2][j][0])][int(clusterTreeEntries['cubes'][iEntry2][j][1])] = 1
            z[j] = clusterTreeEntries['cubes'][iEntry2][j][2]
        if isMuon:
            Zspace = np.linspace(min(z)-5, max(z)+5, 1000)
            jEntry1 = clusterTreeEntries['MuonTreePos'][iEntry2]
            X = muonTreeEntries['mx'][jEntry1] * Zspace + muonTreeEntries['cx'][jEntry1] + 0.5
            Y = muonTreeEntries['my'][jEntry1] * Zspace + muonTreeEntries['cy'][jEntry1] + 0.5

            ax.plot3D(Zspace, X, Y, linewidth=0.7, color='#FF0000', label='{}'.format(eventTwoName))

        ax.voxels(np.array(matrix2), facecolors='#FF000075', edgecolor='#DF010175')

        for m, iC in enumerate(clusterTreeEntries['chanID'][iEntry2]):
            x, y, z, horNotVer = chanIDTranslation(iC)
            if horNotVer and (x == 16):
                Xplusc = np.linspace(-1, 16, 10)
                Yplusc = np.repeat(y + 1, 10)
                Zplusc = np.repeat(z + 1, 10)

                # projection
                rectangle = plt.Rectangle((z, y), 1, 1, fc=colorsElectron[m], ec=colorsElectron[m], alpha=alpha)
                ax.add_patch(rectangle)
                art3d.pathpatch_2d_to_3d(rectangle, z=16, zdir="y")
            elif horNotVer and (x == -1):
                Xplusc = np.linspace(-1, 16, 10)
                Yplusc = np.repeat(y, 10)
                Zplusc = np.repeat(z, 10)

                # projection
                rectangle = plt.Rectangle((z, y), 1, 1, fc=colorsElectron[m], ec=colorsElectron[m], alpha=alpha)
                ax.add_patch(rectangle)
                art3d.pathpatch_2d_to_3d(rectangle, z=-1, zdir="y")

            elif (not horNotVer) and (y == 16):
                Xplusc = np.repeat(x + 1, 10)
                Yplusc = np.linspace(-1, 16, 10)
                Zplusc = np.repeat(z + 1, 10)

                # projection
                rectangle = plt.Rectangle((z, x), 1, 1, fc=colorsElectron[m], ec=colorsElectron[m], alpha=alpha)
                ax.add_patch(rectangle)
                art3d.pathpatch_2d_to_3d(rectangle, z=16, zdir="z")

            else:
                Xplusc = np.repeat(x, 10)
                Yplusc = np.linspace(-1, 16, 10)
                Zplusc = np.repeat(z, 10)

                # projection
                rectangle = plt.Rectangle((z, x), 1, 1, fc=colorsElectron[m], ec=colorsElectron[m], alpha=alpha)
                ax.add_patch(rectangle)
                art3d.pathpatch_2d_to_3d(rectangle, z=-1, zdir="z")

            ax.plot3D(Zplusc, Xplusc, Yplusc, linewidth=0.13, color='#FF0000')  # channels

    ax.set_xbound(0, 50)
    ax.set_ybound(0, 16)
    ax.set_zbound(0, 16)
    ax.set_xlabel('Z (Plane)')
    ax.set_ylabel('X (Cube)')
    ax.set_zlabel('Y (Cube)')
    ax.legend(loc="upper left", bbox_to_anchor=(0.65, 1.05))
    ax.view_init(30, -130)
    ax.grid()
    # add colorbar to figure
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, spacing='uniform', orientation='horizontal',
                                   extend='neither', ticks=[0, 0.25, 0.5, 0.75, 1.0])
    ax1.set_title('Energy [MeV]')
    ax1.set_xticklabels(['{:.0f}'.format(0), '{:.0f}'.format(4000 / 31.5 / 22), '{:.0f}'.format(8000 / 31.5 / 22),
                         '{:.0f}'.format(12000 / 31.5 / 22), '{:.0f}'.format(16000 / 31.5 / 22)])
    ax1.set_position((0.3, 0.1, 0.45, 0.05))

    if numberOfEventsTogether == 1:
        fig.suptitle(r'Visualisation of {} event'.format(eventOneName))
    elif numberOfEventsTogether == 2:
        fig.suptitle(r'Visualisation of {} event with {}'.format(eventOneName, eventTwoName))
    plt.savefig(figName)
    plt.show()

#-------------------------------------------------------------------------------------------		
VisualFromFile(filename='michelEvents.txt', eventNumber=1, isMuon=True, eventOneName='muon', eventTwoName='possible michel electron candidate', figName='muonVisual.png')
VisualFromFile(filename='IBDevents.txt', eventNumber=3, eventOneName='neutron', eventTwoName='positron', figName='IBDevent.png')
