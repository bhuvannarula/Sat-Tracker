'''
Project Sat-Tracker

@team Celestial Observers
@author Bhuvan Narula

Team Members:
    Bhuvan Narula
    Ayush
    Suhana
First Year, IIT Mandi

Author Notes:
- lat, long used in comments is short notation for latitude, longitude
- Program may take some time to load (wait upto 1-2 min)
- When clicking on buttons, windows will close and reopen. DO NOT PANIC or think that 
Program has quit/stopped. Give it a minute.
'''

plot_update_interval = 1000 # in milli sec, edit if required

# Astronomy related Imports
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, TEME, CartesianDifferential, CartesianRepresentation, ITRS
from astropy import units as u

# System Time related Imports
from datetime import datetime, timedelta

# GUI related Imports
from tkinter import Tk, Frame, Label, Button, Entry, StringVar, TOP, BOTH
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk, FigureCanvasTk

# Plotting related Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.basemap import Basemap
from itertools import chain
from matplotlib.animation import FuncAnimation

# Make root GUI window
rootw = Tk()

# Read the data of first satellite in satellite.txt
with open('satellite.txt', 'r', newline = '') as fin:
    t1 = fin.readlines()
    n = t1[0].strip()
    s = t1[1].strip()
    t = t1[2].strip()

# Plotting related variables
marker_color = 'darkslategray'
line_color = 'tab:blue'

# Import satellite data in Satrec
satellite = Satrec.twoline2rv(s, t)

def main_welcome():
    '''
    Tkinter Home Page
    '''
    root = Frame(rootw)
    def gotograph():
        rootw.destroy()
        main_graph()
    def gotofindif():
        root.pack_forget()
        main_where()
    Label(root, text="Welcome to Sat-Track\n", font=('Arial', 25)).pack(side=TOP)
    Label(root, text='We are tracking {} Satellite today.'.format(n)).pack(side=TOP)
    Label(root, text='What do you want to do?\n').pack(side=TOP)
    Button(root, text='Track Satellite', command= gotograph).pack(side=TOP)
    Button(root, text='See if it will pass you!', command= gotofindif).pack(side=TOP)
    root.pack()

def find_if_sat(templat, templon):
    tempskord = SkyCoord(EarthLocation.from_geodetic(templon, templat, 0).itrs)
    def get_speed(t1 = None):
        if not t1:
            t = Time(datetime.now(), format='datetime')
        else:
            t = Time(t1, format='datetime')
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2) # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])

        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        speed1 = ((teme_v.d_x)**2 + (teme_v.d_y)**2 + (teme_v.d_z)**2)**0.5
        return speed1

    def get_skycoord(t1 = None, alt=False):
        if not t1:
            t = Time(datetime.now(), format='datetime')
        else:
            t = Time(t1, format='datetime')
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2) # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])

        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
        itrs = teme.transform_to(ITRS(obstime=t))
        location = itrs.earth_location
        geodat = location.geodetic
        lat0, lon0 = Angle(geodat.lat).degree, Angle(geodat.lon).degree
        skcord = SkyCoord(EarthLocation.from_geodetic(geodat.lon, geodat.lat, geodat.height).itrs) 
        if not alt:
            return skcord
        else:
            return skcord, geodat.height


    def minima():
        ti = datetime.now()
        dt = 30
        dtt = timedelta(seconds=dt)
        cpt = get_skycoord(ti)
        err_dist = get_speed(ti)*(dt*u.s)
        found = False
        pf = []
        while True:
            ti = ti + dtt
            tpt, sat_alt = get_skycoord(ti, alt=True)
            #print(tpt.separation_3d(cpt), err_dist*(0.7))
            tsp = tpt.separation_3d(cpt)
            tsp2 = tpt.separation_3d(tempskord)
            if tsp < err_dist*0.9:
                return False
            elif tsp2 < (sat_alt + 100*u.km):
                pf.append([tsp2, ti, tpt])
                found = True
            else:
                if found:
                    return min(pf, key = lambda var : var[0])
    return minima()


def main_where():
    root = Frame(rootw)
    Label(root, text="Enter your coordinate (Latitude, Longitude) (eg. '0, 0' without quotes): ").pack(side=TOP)
    eb = StringVar(root, value = '')
    Entry(root, textvariable= eb, width=20).pack(side=TOP)
    bt1 = Button(root, text='Submit')
    reslbl = StringVar(root, value = '')
    lbl = Label(root, textvariable=reslbl)
    def cmd1():
        lat, lon = eb.get().split(',')
        lat = float(lat.strip())
        lon = float(lon.strip())
        kk = find_if_sat(lat, lon)
        if not kk:
            reslbl.set('Sorry, Satellite will not be visible at your Location!')
        else:
            reslbl.set('Satellite will be visible, and pass from near you approximately on {}\n(gap < 500 km measured from land to satellite in sky)'.format(kk[1]))
    bt1.configure(command= cmd1)
    bt1.pack(side=TOP)
    lbl.pack(side=TOP)
    root.pack()


def main_graph():
    '''
    Creates Window for Plot of tracking Satellite
    '''
    root = Tk()
    def get_cord(t1 = None, alt = False):
        '''
        Returns tuple containing (latitude, longitude) of Satellite at time 't1'
        If t1 is not specified, System Time is taken
        '''
        if not t1:
            t = Time(datetime.now(), format='datetime')
        else:
            t = Time(t1, format='datetime')
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2) # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])

        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
        itrs = teme.transform_to(ITRS(obstime=t))
        location = itrs.earth_location
        geodat = location.geodetic
        lat0, lon0 = Angle(geodat.lat).degree, Angle(geodat.lon).degree
        if not alt:
            return lat0, lon0
        else:
            return lat0, lon0, geodat.height
    # Create the figure where coordinates will be plotted.
    fig = Figure(figsize=(8, 6), edgecolor='w')
    ax = fig.add_subplot(111)
    # We choose a cylindrical projection of earth
    m = Basemap(projection='cyl', resolution=None,
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, ax=ax)

    def draw_map(m, scale=0.2):
        '''
        Responsible for style of World Map
        '''
        # draw a shaded-relief image
        m.shadedrelief(scale=scale)
        
        # lats and longs are returned as a dictionary
        lats = m.drawparallels(np.linspace(-90, 90, 13))
        lons = m.drawmeridians(np.linspace(-180, 180, 13))

        lat_lines = chain(*(tup[1][0] for tup in lats.items()))
        lon_lines = chain(*(tup[1][0] for tup in lons.items()))
        all_lines = chain(lat_lines, lon_lines)
        
        # setting the style
        for line in all_lines:
            line.set(linestyle='-', alpha=0.3, color='w')

    # Designs the World Map
    draw_map(m)

    def make_path(m, t_show = 120):
        '''
        Responsible for plotting the approximate path which satellite will travel on for next two hours
        t_show : int
        it is minutes after lat, long of satellite should be fetched
        '''
        latt, lonn = [], []
        t1 = datetime.now()
        dt = timedelta(minutes= 1)
        lt1, lg1 = get_cord(t1)
        latt.append(lt1)
        lonn.append(lg1)
        i = 0
        while i < t_show:
            t1 = t1 + dt
            lt1, lg1 = get_cord(t1)
            if lg1*lonn[-1] < -10:
                m.plot(lonn, latt, 'o-', markersize=0, linewidth=1, color=line_color)
                lonn = []
                latt = []
            latt.append(lt1)
            lonn.append(lg1)
            i += 1
        m.plot(lonn, latt, 'o-', markersize=0, linewidth=1, color=line_color)

    # Draw the path
    make_path(m)

    # Initialise the marker which will represent the satellite (a solid circle)
    sctt = m.scatter([], [], latlon = True, marker = 'o', color=marker_color)

    # current lat, long, altitude variable for tkinter
    curdat = StringVar(root, value='')

    def upd_cord(i):
        '''
        Function called to refresh the plot
        '''
        # Fetching the current lat, long of satellite
        lat0, lon0, alt = get_cord(alt = True)
        # Placing the marker at that lat, long
        sctt.set_offsets([lon0, lat0])
        curdat.set('{} Satellite\nLatitude {}\nLongitude {}\nAltitude {}\nTime {}'.format(n, lat0, lon0, alt, datetime.now()))

    # Initialise Canvas for tkinter GUI
    plotcanvas = FigureCanvasTkAgg(fig, root)
    plotcanvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    ani = FuncAnimation(fig, upd_cord, interval = plot_update_interval, blit = False)
    # Initialise Canvas Toolbar for tkinter Canvas
    toolbar = NavigationToolbar2Tk(plotcanvas, root)
    toolbar.update()
    Label(root, textvariable = curdat).pack(side=TOP)
    plotcanvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    plt.show()
    root.pack()
    root.mainloop()

main_welcome()
rootw.mainloop()