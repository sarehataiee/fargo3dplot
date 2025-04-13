import numpy as np
pi= np.pi

class Grid():
    """
    WHAT:  a class that contains the data about the grid
    ========
    INPUT: 
        directory (string): where the output files are located
    ========
    OUTPUTS:
        directory (string): same as directory (just in case!)
        x (array): the azimuthal grid in polar
        y (array): the radial grid w/o ghost cells in polar
        z (array): the vertical grid
        ywg (array): the radial grid w ghost cells
        nx (int): NX
        ny (int): NY
        nz (int): NZ
        xmid (array): center of x
        ymid (array): center of y
        zmid (array): center of z
        nghx(int): number of ghost cells in x
        nghy(int): number of ghost cells in y
        nghz(int): number of ghost cells in z
    """
    def __init__(self, directory):
        # check if the directory is in the proper form
        if directory != './':
            if directory[-1] != '/':
                directory += '/'
        # read x domain
        try:
            x = np.loadtxt(directory+"domain_x.dat")
        except IOError:
            print("IOError with domain_x.dat")
        # read y domain
        try:
            y = np.loadtxt(directory+"domain_y.dat")
        except IOError:
            print("IOError with domain_y.dat")
        # read z domain
        try:
            z = np.loadtxt(directory+"domain_z.dat")
        except IOError:
            print("IOError with domain_z.dat")
        # find the number of ghost cells
        nghx = 0; nghy = 0; nghz = 0
        try:
            data = open(directory+"summary0.dat",'r').readlines()
            for line in data:
                if "NGHY" in line:
                    members = line.split()
                    for each in members:
                        if "NGHY" in each:
                            nghy = int(each[-1])
                if "NGHX" in line:
                    members = line.split()
                    for each in members:
                        if "NGHX" in each:
                            nghx = int(each[-1])
                if "NGHZ" in line:
                    members = line.split()
                    for each in members:
                        if "NGHZ" in each:
                            nghz = int(each[-1])
        except IOError:
            print("IOError with summary0.dat")
        self.directory = directory
        self.x = x
        self.y = y[nghy:-nghy]
        self.z = z[nghz:-nghz]
        self.xwg = x
        self.ywg = y
        self.zwg = z
        self.nghx = nghx
        self.nghy = nghy
        self.nghz = nghz
        self.xmid = 0.5*(self.x[:-1]+self.x[1:])
        self.ymid = 0.5*(self.y[:-1]+self.y[1:])
        self.zmid = 0.5*(self.z[:-1]+self.z[1:])
        if self.xmid.size == 1: 
            self.nx = 0
        else: 
            self.nx = self.xmid.size
        if self.ymid.size == 1: 
            self.ny = 0
        else: 
            self.ny = self.ymid.size
        if self.zmid.size == 1: 
            self.nz = 0
        else: 
            self.nz = self.zmid.size

class ReadField():
    """
    WHAT:  a class contains a given scaler field (polar)
    ========
    INPUTS: 
        directory (string): where the output files are located
        name (string): the field which is going to be read including
                       "dens" for surface density
                       "energy" for energy
                       "vx" 
                       "vy" 
                       "vz" 
                       "press" for pressure
                       
        nout (int): the output number
    ========
    OUTPUT:
        field (array): an (NX,NY,NZ)-shaped array filled with the 2D field
        name (string): name of the field
        nout (int): the output number
        time (float): time in code unit
    """
    def __init__(self, directory, name, nout):
        # check if the directory is in the proper form
        if directory != './':
            if directory[-1] != '/':
                directory += '/'
        # read the field
        self.name = name
        self.nout = nout
        grid = Grid(directory)
        grid_size = np.array([grid.nz, grid.ny, grid.nx])
        non_zeros = grid_size[grid_size.nonzero()]
        if self.name in ["dens", "energy", "vx", "vy", "vz"]:
            try:
                filename = "{}gas{}{}.dat".format(directory,name,nout)
                self.field = np.fromfile(filename).reshape(non_zeros)
            except IOError:
                print("IOError with {}".format(filename))
        if self.name == "press":
            # find adiabatic index
            try:
                data = open(directory+f"summary{nout:d}.dat",'r').readlines()
                for line in data:
                    if "GAMMA" in line:
                        members = line.split()
                        gamma = float(members[-1])
            except IOError:
                print(f"IOError with summary{nout:d}.dat file")
            
            energy = ReadField(directory=directory, name='energy', nout=nout)
            self.field = energy.field * (gamma-1) 
            
        # find the time 
        try:
            data = open(directory+f"summary{nout:d}.dat",'r').readlines()
            for line in data:
                if "simulation time" in line:
                    members = line.split()
                    self.time = float(members[5])
        except IOError:
            print(f"IOError with summary{nout:d}.dat file")
        #"vorticity" for vorticity (only for 2d polar)
        #"vortensity" for vorticity divided by density (only for 2d polar)
        # elif self.name == "vorticity":
        #     # vorticity (= curl of v) is calculated at the lower right corner of each grid
        #     vtheta = ReadField(directory, 'vx', nout, nx, ny).field
        #     vr = ReadField(directory, 'vy', nout, nx, ny).field
        #     grid = Grid(directory)
        #     r = grid.y
        #     theta = grid.x
        #     rmid = grid.ymid
        #     tmid = grid.xmid
        #     dtheta = tmid[1]-tmid[0]  #evenly spaced theta
        #     t2d, r2d = np.meshgrid(theta,r)
        #     t2dm, r2dm = np.meshgrid(tmid,rmid)
        #     # d(rvtheta)/dr
        #     drdvtheta = (r2dm[1:,:]*vtheta[1:,:]-r2dm[:-1,:]*vtheta[:-1,:])/(r2dm[1:,:]-r2dm[:-1,:])
        #     # adding the first ring equals to the second ring
        #     drdvtheta = np.concatenate(([drdvtheta[0,:]],drdvtheta), axis=0)
        #     # d(vr)/dtheta
        #     dthetadvr = (vr[:,1:]-vr[:,:-1])/dtheta
        #     # calculating the last column
        #     dthetadvr = np.concatenate((dthetadvr,np.array([((vr[:,0]-vr[:,-1])/dtheta)]).T), axis=1)
        #     #here is the vorticity
        #     vorticity = (drdvtheta - dthetadvr)/r2dm
        #     # vorticity is calculated now at the grid corners, 
        #     # we need to bring it to the grid centre
        #     # So we do a bilinear interpolation ...
        #     vorticity_c = vorticity[:-1,:-1] * r2d[:-2,:-2] * (t2d[:-2,1:-1]-t2dm[:-1,:-1])  * (r2d[1:-1,:-2]-r2dm[:-1,:-1])
        #     vorticity_c += vorticity[:-1,1:] * r2d[:-2,:-2] * (t2dm[:-1,:-1]-t2d[:-2,:-2])   * (r2d[1:-1,:-2]-r2dm[:-1,:-1])
        #     vorticity_c += vorticity[1:,:-1] * r2d[:-2,:-2] * (t2d[:-2,1:-1]-t2dm[:-1,:-1])  * (r2dm[:-1,:-1]-r2d[:-2,:-2])
        #     vorticity_c += vorticity[1:,1:]  * r2d[:-2,:-2] * (t2dm[:-1,:-1]-t2d[:-2,:-2])   * (r2dm[:-1,:-1]-r2d[:-2,:-2])
        #     # we need to fill the last column
        #     vorticity_c = np.concatenate((vorticity_c,np.array([vorticity_c[:,0]*0]).T), axis=1)
        #     vorticity_c[:,-1] = vorticity[:-1,-1] * r[:-2] * (theta[-1]-tmid[-1])  * (r[1:-1]-rmid[:-1])
        #     vorticity_c[:,-1] += vorticity[:-1,0] * r[:-2] * (tmid[-1]-theta[-2])   * (r[1:-1]-rmid[:-1])
        #     vorticity_c[:,-1] += vorticity[1:,-1] * r[:-2] * (theta[-1]-tmid[-1])  * (rmid[:-1]-r[:-2])
        #     vorticity_c[:,-1] += vorticity[1:,0]  * r[:-2] * (tmid[-1]-theta[-2])   * (rmid[:-1]-r[:-2])
        #     # and low the last row
        #     vorticity_c = np.concatenate((vorticity_c,[vorticity_c[-1,:]]), axis=0)
        #     #
        #     vorticity_c /= r2d[:-1,:-1] * (t2d[:-1,1:]-t2d[:-1,:-1]) * (r2d[1:,:-1]-r2d[:-1,:-1])
        #     self.field = vorticity_c
        # elif self.name == "vortensity":
        #     dens = ReadField(directory, 'dens', nout, nx, ny).field 
        #     vorticity = ReadField(directory, 'vorticity', nout, nx, ny).field
        #     vortensity = vorticity/dens
        #     self.field = vortensity

class ReadPlanet():
    """
    WHAT:  a class gives a planet's quantities. It DOES need the legacy files
    ========
    INPUTS: 
        directory (string): where the output files are located
        nplanet (int): the index of the planet 
    ========
    OUTPUT:
        directory (string): just in case!
        index (int): which of the planets it is
        time (array): time in orbit at 1 unit, read from the orbit file
        semi (array): the planet's semi-major axis
    ========
    METHODS:
        GiveEcc, GiveOrbitalElements, GiveTorquePower, GiveXY, GiveTorqueDensity
    """
    def __init__(self, directory, nplanet):
        # check if the directory is in the proper form
        if directory != './':
            if directory[-1] != '/':
                directory += '/'
        # set the planet's name
        self.directory = directory
        self.index = nplanet
        filename = "{}orbit{}.dat".format(directory,self.index)
        # read the time and semi-major axis from the orbit file
        try:
            time, semi = np.loadtxt(filename, usecols=[0,2], unpack=True)
            self.time = time/2/pi
            self.semi = semi
        except IOError:
            print("IOError with {}".format(filename))
    #
    def GiveEcc(self):
        """
        Returns
        -------
            e: array
               eccentricity
        """
        filename = "{}orbit{}.dat".format(self.directory,self.index)
        e = np.loadtxt(filename, usecols=[1])
        return e
    #
    def GiveOrbitalElements(self):
        """
        Returns
        -------
        M : array
            mean anomaly
        paromega : array
            argument of pericentre
        lambda_p : array
            mean longitude that is M+paromega (-pi <lambda_p< pi)
        """
        filename = "{}orbit{}.dat".format(self.directory,self.index)
        M, paromega = np.loadtxt(filename, usecols=[3,5], unpack=True)
        lambda_p = M + paromega
        while (lambda_p.min() < -pi):
            np.putmask(lambda_p, lambda_p < -pi*np.ones(lambda_p.size), lambda_p + 2.*pi)
        while (lambda_p.max() > pi):
            np.putmask(lambda_p, lambda_p > pi*np.ones(lambda_p.size), lambda_p - 2.*pi)
        return M, paromega, lambda_p
    #
    def GiveTorquePower(self):
        """
        Read the tqwk file that is one of the legacy files. Thus, the line
            FARGO_OPT += -DLEGACY
        must be in your .opt file

        Returns
        -------
        total_torque: array
            torque from the inner and outer disc 
        total_power: array
            power from the inner and outer disc 
        torq_in : array
            torque from the inner disc 
        torq_out : array
            torque from the outer disc 
        pow_in : array
            power from the inner disc 
        pow_out : array
            power from the outer disc 
        """
        filename = "{}tqwk{}.dat".format(self.directory,self.index)
        torq_in, torq_out, pow_in, pow_out = np.loadtxt(filename, usecols=[1,2,5,6], unpack=True)
        total_torque = torq_in+torq_out
        total_power = pow_in+pow_out
        return total_torque, total_power, torq_in, torq_out, pow_in, pow_out
    #
    def GiveTorqueOnly(self):
        """
        Read the torq_planet file that is one of the monitoring files. Thus the line
            MONITOR_Y  = TORQ
        must be in your .opt file.

        Returns
        -------
        date: array
            the time in the unit of 2*pi
        total_torque: array
            total torque on the planet
        """
        filename = "{}monitor/gas/torq_planet_{}.dat".format(self.directory,self.index)
        date, total_torque = np.loadtxt(filename, unpack=True)
        return date, total_torque
    #
    def GiveXY(self):
        """
        Read the position of the planet from the bigplanet file that is one of the legacy files. 
        Thus, the line
            FARGO_OPT += -DLEGACY
        must be in your .opt file

        Returns
        -------
        x : array
            azimuth of the planet
        y : array
            radial position of the planet

        """
        filename = "{}bigplanet{}.dat".format(self.directory,self.index)
        x, y = np.loadtxt(filename, usecols=[1,2], unpack=True)
        return x,y
    #
    def GiveTorqueDensity(self, what_time, ny):
        """
        Read the torq_1d_Y_raw_planet file that is a monitiring file. Thus the line
            MONITOR_Y_RAW  = TORQ
        must be in your .opt file.

        Parameters
        ----------
        what_time : float
            the time you want to have the torque density (in orbit)
        ny : int
            NY

        Returns
        -------
        torque_dens: a 1D NY-size array
            the torque from each radius on the planet

        """
        filename = "{}monitor/gas/torq_1d_Y_raw_planet_{}.dat".format(self.directory,self.index)
        data = np.fromfile(filename)
        try:
            i_time = np.where(self.time >= what_time)[0][0]
        except IndexError:
            print("Oops! It seems your time has not reached, but you can have the first output :-)")
            i_time = 0
        torque_dens = data[ny*i_time:ny*(i_time+1)]
        return torque_dens
