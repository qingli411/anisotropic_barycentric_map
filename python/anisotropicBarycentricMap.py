import numpy as np
import matplotlib.pyplot as plt

class AnisotropyTensor(object):
    """ ReynoldsStress object
    """

    def __init__(self, uu, vv, ww, uv, uw, vw, ):
        """Initialize an anisotropy tensor object from the six components of Reynolds stress.

        """
        data = {}
        data['uu'] = np.array(uu)
        data['vv'] = np.array(vv)
        data['ww'] = np.array(ww)
        data['uv'] = np.array(uv)
        data['uw'] = np.array(uw)
        data['vw'] = np.array(vw)
        self.size = data['uu'].size
        for key in data.keys():
            assert data[key].ndim <= 1, "Each component of the Reynolds stress should be a scalar or a 1D array."
            assert data[key].size == self.size, "All components of the Reynolds stress should have the same size."
        tke = 0.5*(data['uu']+data['vv']+data['ww'])
        self.data = np.zeros([self.size, 3, 3])
        self.data[:,0,0] = 0.5*data['uu']/tke - 1/3
        self.data[:,1,1] = 0.5*data['vv']/tke - 1/3
        self.data[:,2,2] = 0.5*data['ww']/tke - 1/3
        self.data[:,0,1] = 0.5*data['uv']/tke
        self.data[:,0,2] = 0.5*data['uw']/tke
        self.data[:,1,2] = 0.5*data['vw']/tke
        self.data[:,1,0] = 0.5*data['uv']/tke
        self.data[:,2,0] = 0.5*data['uw']/tke
        self.data[:,2,1] = 0.5*data['vw']/tke

    def __repr__(self):
        """Formatted print

        """
        summary = ['{}'.format(self.data[0,:,:])]
        if self.size > 1:
            summary.append('...')
            summary.append('{}'.format(self.data[-1,:,:]))
        return '\n'.join(summary)

    def barycentric_coord(self):
        """Get the coordinates in the barycentric map

        See Banerjee et al., 2007

        """
        c = np.zeros([self.size, 3])
        for i in np.arange(self.size):
            a = self.data[i,:,:]
            eigenvalues, _ = np.linalg.eig(a)
            idx = eigenvalues.argsort()[::-1]
            l = eigenvalues[idx]
            c[i,0] = l[0]-l[1]
            c[i,1] = 2*(l[1]-l[2])
            c[i,2] = 3*l[2]+1
        return c

    def vector_direction_coord(self):
        """ Get the coordinates in the vector direction map

        See Li and Fox-Kemper, 2020

        """
        cmin = np.zeros([self.size, 3])
        cmax = np.zeros([self.size, 3])
        for i in np.arange(self.size):
            a = self.data[i,:,:]
            eigenvalues, eigenvectors = np.linalg.eig(a)
            idx = eigenvalues.argsort()[::-1]
            l = eigenvalues[idx]
            v = eigenvectors[:,idx]
            eta1 = l[0]-l[1]
            eta2 = 2*(l[1]-l[2])
            if eta1 >= eta2:
                cmin[i,:] = np.nan
                rv = np.sqrt(v[0,0]**2+v[1,0]**2)
                phi = np.arctan2(v[2,0], rv)
                theta = np.arctan2(v[1,0], v[0,0])
                rr = 1-np.abs(phi/np.pi*2)
                cmax[i,0] = rr*np.cos(theta)*np.sign(phi)
                cmax[i,1] = rr*np.sin(theta)*np.sign(phi)
                cmax[i,2] = eta1
            else:
                cmax[i,:] = np.nan
                rv = np.sqrt(v[0,2]**2+v[1,2]**2)
                phi = np.arctan2(v[2,2], rv)
                theta = np.arctan2(v[1,2], v[0,2])
                rr = 1-np.abs(phi/np.pi*2)
                cmin[i,0] = rr*np.cos(theta)*np.sign(phi)
                cmin[i,1] = rr*np.sin(theta)*np.sign(phi)
                cmin[i,2] = eta2
        return [cmax, cmin]

    def plot_anisotropic_barycentric_map(self, ax=None, **kwargs):
        """Plot anisotropic barycentric map

        """
        if ax is None:
            ax = plt.gca()
        if not ax.has_data():
            plot_anisotropic_barycentric_map_base(ax)
        xc = np.array([0.5, -0.5, 0])
        yc = np.array([0, 0, np.sqrt(3)*0.5])
        c = self.barycentric_coord()
        xx = np.dot(xc, c.transpose())
        yy = np.dot(yc, c.transpose())
        im = ax.scatter(xx, yy, zorder=3, **kwargs)
        return im

    def plot_eigen_vector_direction_max_min(self, ax=None, c=None, cmap='rainbow'):
        """Plot the direction of the eigen vectors corresponding to the maximum and minimum eigen values.

        """
        if ax is None:
            ax = plt.gca()
        if not ax.has_data():
            plot_vector_direction_base(ax)
        cmax, cmin = self.vector_direction_coord()
        bsize = 200
        cmap_colors = plt.cm.get_cmap(cmap)
        frac = c/np.abs(c).max()
        colors = cmap_colors(1+frac)
        im = ax.scatter(cmin[:,0], cmin[:,1], s=cmin[:,2]*bsize, c=c, marker='x', cmap=cmap, zorder=3)
        ax.scatter(cmax[:,0], cmax[:,1], s=cmax[:,2]*bsize, marker='o', edgecolors=colors, facecolors='none', linewidth=1.5, zorder=3)
        ax.scatter(1, 1.1, s=bsize, marker='o', edgecolor='k', facecolor='none', linewidth=1.5)
        ax.scatter(1, 0.9, s=bsize, c='k', marker='x')
        return im

def plot_anisotropic_barycentric_map_base(ax):
    # plot a triangle
    xc = np.array([0.5, -0.5, 0])
    yc = np.array([0, 0, np.sqrt(3)*0.5])
    for i in np.arange(3):
        ip1 = (i+1)%3
        ax.plot([xc[i], xc[ip1]], [yc[i], yc[ip1]], 'k', linewidth=2)
    # add grid
    nsp = 5
    lc = np.abs(xc[1]-xc[0])
    dc = lc/nsp
    cl = np.zeros([3*(nsp-1), 3])
    for i in np.arange(3):
        ip1 = (i+1)%3
        for j in np.arange(nsp-1):
            k = i * (nsp-1) + j
            cl[k,i] = dc * (j+1)
            cl[k,ip1] = dc * (nsp-j-1)
    xl = np.dot(xc, cl.transpose())
    yl = np.dot(yc, cl.transpose())
    nl = xl.size
    for i in np.arange(3):
        ip1 = (i+1)%3
        for j in np.arange(nsp-1):
            k = i * (nsp-1) + j
            kp = (ip1 * (nsp-1) + nsp - j - 2) % nl
            ax.plot([xl[k], xl[kp]], [yl[k], yl[kp]], '--k', linewidth=0.75)
    # plain strain limit
    c1ps = np.array([1/3, 2/3, 0])
    x1ps = np.dot(xc, c1ps.transpose())
    y1ps = np.dot(yc, c1ps.transpose())
    c2ps = np.array([0, 0, 1])
    x2ps = np.dot(xc, c2ps.transpose())
    y2ps = np.dot(yc, c2ps.transpose())
    ax.plot([x1ps, x2ps], [y1ps, y2ps], '-k', linewidth=0.75)
    # add labels
    labels = ['1 comp', '2 comp', '3 comp']
    lbpx = xc
    dshift = 0.05
    lbpy = [yc[0]-dshift*lc, yc[1]-dshift*lc, yc[2]+dshift*lc]
    for i in np.arange(3):
        ax.text(lbpx[i], lbpy[i], labels[i], ha='center', fontsize=12)
    label_side1 = 'Axisymmetric expansion'
    label_side2 = 'Axisymmetric contraction'
    label_side3 = 'Two component'
    label_ps = 'Plane strain'
    ax.text((xc[0]+xc[2])/2, (yc[0]+yc[2])/2+0.08*lc, label_side1, ha='center', va='center', rotation=-60)
    ax.text((xc[1]+xc[2])/2, (yc[1]+yc[2])/2+0.08*lc, label_side2, ha='center', va='center', rotation=60)
    ax.text((xc[0]+xc[1])/2, (yc[0]+yc[1])/2-0.04*lc, label_side3, ha='center', va='center')
    ang = np.arctan2(np.abs(y1ps-y2ps), np.abs(x1ps-x2ps))/np.pi*180
    ax.text((x1ps+x2ps)/2, (y1ps+y2ps)/2+0.16*lc, label_ps, ha='center', va='center', rotation=ang)

    # set axis
    ax.axis('equal')
    ax.axis('off')

def plot_vector_direction_base(ax):
    # plot a circle
    r = 1.
    theta = np.deg2rad(np.linspace(0,360,361))
    xcircle = r*np.cos(theta)
    ycircle = r*np.sin(theta)
    ax.plot(xcircle, ycircle, color='gray', linewidth=2)
    # x-y angles
    nr = 12
    for i in np.arange(nr):
        phi = np.pi*2/nr*i
        xref = r*np.cos(phi)
        yref = r*np.sin(phi)
        ax.plot([0, xref], [0, yref], color='gray', linestyle='--', linewidth=1)
        if i > 0 and i != nr/4:
            di = i/nr*360
            if di > 180:
                di = di-360
            rlabel = '{:3.0f}$^\circ$'.format(di)
            xlb = 1.15*xref
            ylb = 1.15*yref
            ax.text(xlb, ylb, rlabel, color='gray', ha='center', va='center')
    # z angle
    nz = 5
    for i in np.arange(nz):
        ri = r/(nz+1)*(i+1)
        xref = ri*np.cos(theta)
        yref = ri*np.sin(theta)
        ax.plot(xref, yref, color='gray', linestyle='--', linewidth=1)
        zlabel = '{:2.0f}$^\circ$'.format(90*(1-ri/r))
        xlb = -ri-0.1
        ylb = -0.08
        ax.text(xlb, ylb, zlabel, color='gray', fontsize=8)
    # plot x-y axes
    xx = 1.1
    yy = 1.1
    ax.arrow(0, 0, xx, 0, linewidth=2, head_width=0.04, color='k', zorder=3)
    ax.arrow(0, 0, 0, yy, linewidth=2, head_width=0.04, color='k', zorder=3)
    ax.text(xx, -0.1, '$x$', ha='center', va='center', fontsize=12)
    ax.text(-0.1, yy, '$y$', ha='center', va='center', fontsize=12)
    ax.text(0.06, 0.08, '$z$', ha='center', va='center', fontsize=12)

    # set axis
    ax.axis('equal')
    ax.axis('off')
