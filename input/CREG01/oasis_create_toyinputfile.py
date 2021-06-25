#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Valerie Garnier, Joris Pianezze '
__email__ = 'vgarnier@ifremer.fr'
__date__ = '2015-03-12'
__doc__ = """  DESCRIPTION \n   =========== \n

    Creation of grid and restart files for oasis use

    Assumption : the netCDF file correspond to the grid exchanged by OASIS (.i.e first row and column removed for MARS)

    The netCDF files must contains :
        The variables of interest SST SSH SSU SSV...
        The mask or bathymetry
        At least the longitude and latitude at t-location

    REQUIRED :
    ----------

        an [list of] input file
        --input_model : model from which data come from [MARS, CROCO, MNH or WW3
        --output_model : model which will read the data (via OASIS) to start the run [MARS, CROCO, MNH, WW3 or TOY]

    OPTIONAL :
    ----------

        --gridfile : to be specified if the grid is in a separate netCDF file
        --indir : to specify the directory where are the data to be processes
        --outdir : to specify the directory where you want to save the figures
        --add_lonlattime : to specify the real longitude, latitude and time of records
        -p : to print intermediate logs

    DVPER use : the trick is to file properly the classes
        - Model (parent of the other classes and defines all attributes of all the fields)
        - WW3
        - MesoNH
        - Croco
        - Mars
        - Toy
        - Oasis (to use export files of OASIS)
        The key of the dictionaries defined the fields are fixed because allow the recognize the fields whatever the model
        the dvper should change only the value or add new fields.

        Fields to  make sure nothing is going wrong with new changes under /home/caparmor-work/marsdev/NONREGRESSION

    EXAMPLE :
        MARS :
        oasis_create_rstrt.py champs_MENOR_MARS_VB2142.nc --input_model MARS -p --create_grid --output_model TOY

        CROCO :
        oasis_create_rstrt.py roms_his_iroise.nc --indir PATH --gridfile roms_grid_iroise.nc --input_model CROCO -p --create_grid --output_model CROCO

        WW3 :
        oasis_create_rstrt.py ww3_iroise.nc --input_model WW3 -p --create_grid --output_model WW3

        MNH :
        oasis_create_rstrt.py TEST1.1.09-02.024.nc --create_grid --indir PATH --input_model MNH -p --gridfile TEST1.1.09-02.024.nc --output_model TOY --mnhtimeref "20110902 000000" ----mnhtimerec 3600
        oasis_create_rstrt.py TEST1.1.09-02.024sf2.nc --create_grid --indir PATH --input_model MNH -p --gridfile TEST1.1.09-02.024.nc --output_model TOY --mnhtimeref "20110902 000000" ----mnhtimerec 3600
        connect_files.py rstrt_from_MNH_to_TOY_fromfile_TEST1.1.09-02.024.nc,rstrt_from_MNH_to_TOY_fromfile_TEST1.1.09-02.024sf2.nc,rstrt_from_MNH_to_TOY_fromfile_TEST1.1.09-02.024sf2.nc --var TOY_FMSU,TOY_FMSV --outfile rstrt_from_MNH_to_TOY_fromfile_allTEST1.1.09-02.024.nc

        OASIS  (on suppose la grille connue) :
        oasis_create_rstrt.py oasis.nc --outdir PATH --input_model OASIS -p --gridfile TEST1.1.09-02.024.nc --output_model TOY


    Commentaire pour futur : il me reste à faire

        - sauver les temps (les lire quand ils existent et sont superieur à 1 (ou ont time_origine), les recuperer a partir du nom de fichier pour MNH
        - separer partie grid et rstrt ?

    Bilan
        - MARS / MARS et MARS / TOY nickel. Encore heureux vu que le codage est parti de mars...
        - CROCO / CROCO et CROCO / TOY a priori OK.
        - MNH / MNH et MNH / TOY a priori OK Utiliser connect_files.py pour concatener les champs venant de sf2
        - WW3 / WW3 et WW3 / TOY a priori OK
             * hum sans doute rajouter une routine pour passer de FP a TP, recuperer la norme de ubr, voir les cos sin de dir aussi. Bizarre, je ne devrais pas avoir a faire cela. A voir avec la pratique.
        - OASIS / TOY Utiliser connect_files.py pour concatener les champs sortant de OASIS

   END DESCRIPTION \n   ===============
   """


import os, time, sys, datetime
import argparse
import netCDF4
import numpy as np
from itertools import chain
import toolsoasis


def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

def set_wetdry_mask(filename,typecode):

    """ Get the wet-dry mask : 0 for land, 1 for sea

    :param filename: contains the data
    :param typecode: class describing the possible content of the file (according to MASR or CROCO or...)
    :return:a dictionnary with a few of the following keys : drywet,drywet_x,drywet_y
    """

    wetdry = {}

    try:

        if typecode.model == "MARS":
            for im,name in enumerate(typecode.topo):

                # Read the bathymetry
                bathy = input.get_topo(filename,name,typecode)

                # Check and set where is dry
                mask = np.ones(bathy.shape)
                mask = np.where(bathy[:] > -50., mask, 0.) 

                # Add the mask
                if name.endswith("0"): wetdry.update( { "drywet" : mask})
                elif name.lower().endswith("x"): wetdry.update( { "drywet_x" : mask})
                elif name.lower().endswith("y"): wetdry.update( { "drywet_y" : mask})

        if typecode.model == "CROCO":
            for im,name in enumerate(typecode.mask):

                # Read the bathymetry
                mask = input.get_topo(filename,name,typecode)

                # Add the mask
                if name.endswith("o"): wetdry.update( { "drywet" : mask})
                elif name.lower().endswith("u"): wetdry.update( { "drywet_x" : mask})
                elif name.lower().endswith("v"): wetdry.update( { "drywet_y" : mask})

        if typecode.model == "WW3":
            for im,name in enumerate(typecode.mask):

                # Read the mask
                bathy = input.get_topo(filename,name,typecode)

                # Check and set where is dry
                mask = np.ones(bathy.shape)
                mask = np.where(bathy[:] != 8., mask, 0.) 

                # Add the mask
                wetdry.update( { "drywet" : mask})

        if typecode.model == "MNH" or typecode.model == "OASIS":
            for im,name in enumerate(typecode.topo):

                # Read the bathymetry
                bathy = input.get_topo(filename,name,typecode)

                # Check and set where is dry
                mask = np.ones(bathy.shape)
                mask = np.where(bathy[:] == 0., mask, 0.) 

                # Add the mask
                wetdry.update( { "drywet" : mask})

    except Exception, message:
        print "    The %s %s field is not available. None evaluation of the wet-dry mask at this location" %(typecode.model,name)
        print "    Error due to %s" %message

    return wetdry


if __name__ == '__main__' :

    # Parser to specify the choice of the user through arguments given to the script
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog = 'oasis_create_rstrt.py',
        description=__doc__,
        usage="%(prog)s infile.nc --input_model MARS --output_model TOY",
        version=__date__
        )

    parser.add_argument("infile", help="NetCDF files from where you will extract the field",default=None)
    parser.add_argument("--gridfile", help="NetCDF grid file read by OASIS (for the spatial interpolation)",default=None)

    parser.add_argument("--input_model", help="model whose data come from [MARS, CROCO, MNH, WW3 or OASIS]",default=None)
    parser.add_argument("--output_model", help="model which will read the data (via OASIS) to start the run [MARS, CROCO, MNH, WW3 or TOY]",default=None)

    parser.add_argument("--create_grid", action='store_true',help="ask for the creation of the OASIS grid")
    parser.add_argument("--add_lonlattime", action='store_true',help="ask for the specitication of longitude, latitude and time records")
    parser.add_argument('--indir', help='Directory where are the data', default="./")
    parser.add_argument('--outdir', help='Directory where are saved the figures', default="./")

    parser.add_argument("--mnhtimeref", help="time origin of the netCDF file (MesoNH purpose) [default: %(default)s]",default="20110902 000000")
    parser.add_argument("--mnhtimerec", help="time [s] corresponding to the MesoNH netCDF file [default: %(default)s]",default=0)

    parser.add_argument("-p","--print", dest='debug', action='store_true',help="print intermediate variables")


    options = parser.parse_args()
    if options.debug: print "\n options summary %s \n" %options


    # checking of the options
    # !!!!!!!!!!!!!!!!!!!!!!!

    input_model = options.input_model
    output_model = options.output_model
    if options.input_model not in ["MARS","CROCO","MNH","WW3","OASIS"]: parser.error('Specify --input_model [MARS, CROCO, MNH, WW3 or OASIS]')
    if options.output_model not in ["MARS","CROCO","MNH","WW3","TOY"]: parser.error('Specify --output_model [MARS, CROCO, MNH, WW3 or TOY]')

    if not options.gridfile: gridfile = os.path.join(options.indir,options.infile.split(',')[0])
    else: gridfile = os.path.join(options.indir,options.gridfile)

    if options.input_model == "OASIS" and options.create_grid:
        parser.error('with --input_model OASIS you cannot estimate the grid (which already exist anyway)')

    # Save the command line
    # !!!!!!!!!!!!!!!!!!!!!

    f = open(".".join((".".join(('cmdline',__file__.split('/')[-1].replace('.py',''))),'log')) , "a" )
    f.write( time.asctime() +'\n' )
    f.write( __file__ +' with the following options:\n' )
    f.write(str(options)+'\n')
    f.write('\n'*2)
    f.close()

    # Processing
    # !!!!!!!!!!

    # Initialization
    if options.input_model == "MARS": input = toolsoasis.Mars()
    elif options.input_model == "CROCO": input = toolsoasis.Croco()
    elif options.input_model == "MNH": input = toolsoasis.MesoNH()
    elif options.input_model == "WW3": input = toolsoasis.WW3()
    elif options.input_model == "OASIS": input = toolsoasis.Oasis()

    if options.output_model == "MARS": output = toolsoasis.Mars()
    elif options.output_model == "CROCO": output = toolsoasis.Croco()
    elif options.output_model == "MNH": output = toolsoasis.MesoNH()
    elif options.output_model == "WW3": output = toolsoasis.WW3()
    elif options.output_model == "TOY": output = toolsoasis.Toy()


    # Grid
    # ----

    print "\n"
    print "Read the grid file %s" %gridfile

    # Read the grid fields
    gridnc = netCDF4.Dataset(gridfile)

    # longitude and latitude
    lon = input.get_lonlat(gridnc,input.longitude)
    lat = input.get_lonlat(gridnc,input.latitude)
    [lon,lat] = input.check_lonlat_dimension(lon,lat,input)

    # Close the input file
    gridnc.close()


    # Variables
    # ---------

    for ifile,file in enumerate(options.infile.split(',')):

        print "\n"
        print "Read the file %s" %os.path.join(options.indir,file)

        # Read the fields of interest
        infile = netCDF4.Dataset(os.path.join(options.indir,file))

        # get the available dimensions in the file
        axes = []
        try:
            if infile.dimensions[ input.dims["time"]["name"] ]:
                ntime = int( "".join([el for el in str(infile.dimensions[ input.dims["time"]["name"] ])[-10:] if hasNumbers(el) ]) )
                axes.append(('t',ntime))
        except:
            print "    No time axis in the input file %s" %file
        try:
            if infile.dimensions[ input.dims["level"]["name"] ]: pass  #axes.append('z')
        except:
            print "    No level axis in the input file %s" %file

        # exchanged fields (as a whole)
        field = []
        # get the current variable names used in the codes
        name_keys = input.field.keys()
        current_names = []
        current_names.extend([input.field[name_keys[i]]["name"] for i in range(len(name_keys))])
        current_names = list(chain.from_iterable(current_names))

        if options.debug: print " Variables currently available in a %s file : %s" %(input.model,current_names)

        if options.debug: print " Variables available actually %s" %infile.variables.keys()

        # read the field
        list_var = list( set(current_names) & set(infile.variables.keys()) )
        for fieldname in list_var:
            print "    Read the variable %s" %fieldname
            field.append(dict([(fieldname,input.get_var(infile,fieldname,input))]))
#            if input.model=="MARS":   # car pb de masque avec le TOY (en fait ne sert a rien mais modification de fillvalue plus efficace
#                if fieldname in input.field["ssu"]["name"]:
#                    field[-1][fieldname].mask = np.ma.nomask
                    #field[-1][fieldname][:][field[-1][fieldname][:]==input.fillvalue] = 0.
#                    field[-1][fieldname][:][field[-1][fieldname][:] > 90.0] = 0.
#                if fieldname in input.field["ssv"]["name"]:
#                    field[-1][fieldname].mask = np.ma.nomask
#                    field[-1][fieldname][:][field[-1][fieldname][:] > 90.0] = 0.
            if input.model=="MNH":
                if fieldname == "GFLUX_SEA":
                    # Estimate the Net solar heat flux
                    print "    Estimate and add the variable %s" %input.field["swnet"]['name'][0]
                    swu = input.get_var(infile, input.field["swu"]['name'][0] ,input)
                    swd = input.get_var(infile, input.field["swd"]['name'][0] ,input)
                    tmp = swd - swu
                    ctmp = np.ma.masked_where( tmp == 0.0, tmp )
                    fieldname = "swnet"
                    field.append(dict([(fieldname,ctmp)]))
                    # Estimate the Non solar heat flux
                    print "    Estimate and add the variable %s" %input.field["heat"]['name'][0]
                    hftot = input.get_var(infile, input.field["hflux"]['name'][0] ,input)
                    tmp = hftot - (swd-swu)
                    tmp[tmp >= input.fillvalue] = input.fillvalue
                    ctmp = np.ma.masked_where( tmp >= input.fillvalue , tmp )
                    fieldname = "heat"
                    field.append(dict([(fieldname,ctmp)]))
                    # Estimate the Evaporation
                    print "    Estimate and add the variable %s" %input.field["evap"]['name'][0]
                    fieldname = "evap"
                    hflat = input.get_var(infile, input.field["hflat"]['name'][0] ,input)
                    ctmp = np.ma.masked_where( hflat >= input.fillvalue , hflat/2.5e+6 )
                    field.append(dict([(fieldname,ctmp)]))
            if fieldname == "dir":
                tmp = input.get_var(infile,fieldname,input)
                fieldname = "cdir"
                print "    Add the variable %s" %fieldname
                ctmp = np.ma.masked_where(   tmp.mask, np.cos( np.deg2rad(tmp) )  )
                field.append(dict([(fieldname,ctmp)]))
                fieldname = "sdir"
                print "    Add the variable %s" %fieldname
                stmp = np.ma.masked_where(   tmp.mask, np.sin( np.deg2rad(tmp) )  )
                field.append(dict([(fieldname,stmp)]))
                del ctmp,stmp
            if fieldname == "fp":
                tmp = input.get_var(infile,fieldname,input)
                fieldname = "tp"
                print "    Add the variable %s" %fieldname
                ctmp = np.ma.masked_where(   tmp.mask, 1./tmp  )
                field.append(dict([(fieldname,ctmp)]))
        if all((w in list_var for w in ["uubr","vubr"])):
                utmp = input.get_var(infile,"uubr",input)
                vtmp = input.get_var(infile,"vubr",input)
                tmp = (utmp**2 + vtmp**2)**0.5
                field.append(dict([("ubr",tmp)]))
                del tmp,utmp,vtmp

        # Add the wet-dry masks
        wetdry_dict = set_wetdry_mask(infile,input)

        # Save the input file
        outfile_name = "_".join(("input_from_%s_to_%s" %(options.input_model,options.output_model),"fromfile",os.path.basename(file)))
        outfile_name = os.path.join(options.outdir,outfile_name)
        axis = [] # None T oz Z axis
        if len(axes) > 0: axis=axes[0]

        if options.add_lonlattime:
            dt_origin_WW3 = datetime.datetime.strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
            if input.model == "MARS":
                timevar = infile.variables['time']
                timevalue = infile.variables['time'][:]
                dt_origin = datetime.datetime.strptime(timevar.time_origin, "%d-%b-%Y %H:%M:%S")
                set_varaxes = timevalue + int( dt_origin.strftime("%s") ) - int( dt_origin_WW3.strftime("%s") )
            elif input.model == "MNH":
                DELTAT_DATE = float(options.mnhtimerec)
                DT_ORIGIN_SIMU = datetime.datetime.strptime(options.mnhtimeref, "%Y%m%d %H%M%S")
                SECONDS_SINCE_REF = (DT_ORIGIN_SIMU - dt_origin_WW3).total_seconds()
                print "SECONDS_SINCE_REF",SECONDS_SINCE_REF
                set_varaxes = [DELTAT_DATE+SECONDS_SINCE_REF]
                print len(set_varaxes)
            elif input.model == "WW3":
                timevar = infile.variables['time']
                timevalue = infile.variables['time'][:]*86400
                dt_origin = datetime.datetime.strptime(timevar.units[11:], "%Y-%m-%d %H:%M:%S")
                set_varaxes = timevalue + int( dt_origin.strftime("%s") ) - int( dt_origin_WW3.strftime("%s") )
        else: set_varaxes = False
        output.save_inputfile(lon,lat,field,axis,wetdry_dict,outfile_name,input,output,set_varaxes=set_varaxes)

        # Close the input file
        infile.close()


    print "\n"
    print "That's it, guy"
