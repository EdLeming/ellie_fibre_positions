'''
A short script to calculate the TELLIE fibre install positions
using the install spreadsheet (used by the installation team)
and a text file containing the hex panel x,y,z coordinates as
inputs.

The script additionally creates new .json files for where the
tellie fibre tables have been updated with three new fields.

1) psup_panel:     The panel on which a fibre is installed
2) pointing_angle: The pointing angle of the fibre relative to
                   the pointing vector of the panel it's
                   installed on.
3) fibre_status:   A status flag, where: 
                   0 : Installed, OK
                   1 : Installed, low transmittance
                   2 : Installed, broken
                   3 : Not installed, OK
                   4 : Not installed, low transmittance
                   5 : Not installed, broken

The updated files will eventually be copied into TELLIE.ratdb,
stored in the rat/data directory.

Author: Ed Leming
Date  : 31/08/2016
'''

import rat
import ROOT

import optparse
import csv
import numpy as np
import sys 

def read_install_table(fname):
    '''Read in relavent fields from Sofia's install table
    '''
    fibres = []
    nodes = {}
    pmt_hex = {}
    neighbour_hex = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            # Read in the node info - there's one wild one for FT001
            # which isn't a number, just continue past this one.
            try:
                int(row[0])
            except:
                continue
            
            node = int(row[0])
            # Read in fibres, make sure to get the ones wrongly installed!
            if row[3] == '':
                fibre = row[2]
            else:
                fibre = row[3]
            if fibre[:2] == "FT":
                fibres.append("%sA" % fibre)
                fibres.append("%sB" % fibre)
                nodes["%sA" % fibre] = node
                nodes["%sB" % fibre] = node
                # Read in plate positions, make sure to get the ones wrongly installed!
                if row[12] == '':
                    pmt_hex["%sA" % fibre] = row[8]
                    pmt_hex["%sB" % fibre] = row[8]
                    neighbour_hex["%sA" % fibre] = row[9]
                    neighbour_hex["%sB" % fibre] = row[9]
                else: 
                    pmt_hex["%sA" % fibre] = row[12]
                    pmt_hex["%sB" % fibre] = row[12]
                    neighbour_hex["%sA" % fibre] = row[13]
                    neighbour_hex["%sB" % fibre] = row[13]
            else: 
                fibres.append(fibre)
                nodes[fibre] = node
                # Read in plate positions, make sure to get the ones wrongly installed!
                if row[12] == '':
                    pmt_hex[fibre] = row[8]
                    neighbour_hex[fibre] = row[9]
                else: 
                    pmt_hex[fibre] = row[12]
                    neighbour_hex[fibre] = row[13]

    return nodes, fibres, pmt_hex, neighbour_hex

def read_PMT_coordinates(fname):
    '''Get global positions of PMT Hex cells
    '''
    cells = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row == []:
                continue
            row = filter(None, row)
            try:
                cells[row[1]] = [float(row[3]), float(row[4]), float(row[5])]
            except:
                raise
    return cells
    
def get_pmt_coordinates(host_hex, neighbour_hex, fname):
    '''Use the PMT hex cell file to find PMT positions and return
    the relavent three vectors for the plate position calculation
    '''
    # Some brute force for correcting fields read from .csv
    if "none" in neighbour_hex:
        neighbour_hex = host_hex

    panel_pmts = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # Skip header line
        for row in reader:
           # Skip over empty lines
            if row == []:
                continue
            # Remove empty strings
            row = filter(None, row)
            # Get three vectors for host and neighbour cells + record the hex_panel name
            if row[1][:4] == host_hex:
                host_vec = ROOT.TVector3(float(row[3]), float(row[4]), float(row[5]))
                panel_name = row[6][:4]
            if row[1][:4] == neighbour_hex:
                neighbour_vec = ROOT.TVector3(float(row[3]), float(row[4]), float(row[5]))
        # If no neighbour listed then panel is directly above host
        try:
            neighbour_vec
        except:
            neighbour_vec = host_vec

        # Now find all cells in that panel!
        f.seek(0)
        next(reader, None) # Skip header line
        for row in reader:
            # Skip over empty lines
            if row == []:
                continue
            # Remove empty strings
            row = filter(None, row)
            if row[6][:4] == panel_name:
                panel_pmts[int(row[6][-2:])] = [float(row[3]), float(row[4]), float(row[5])]

    # Find central PMT in this hex panel.
    if len(panel_pmts) > 9:
        central_pmt = panel_pmts[10]
    else:
        central_pmt = panel_pmts[4]
    central_vec = ROOT.TVector3(central_pmt[0], central_pmt[1], central_pmt[2])
    return central_vec, host_vec, neighbour_vec, int(panel_name[1:])

def calc_fibre_placement(host_hex, neighbour_hex, fname):
    '''Use the install hex and neighbour hex to calculate the plate install position
    '''
    # Define hex cell & plate parameters
    hex_size = 17.05                             # [cm]
    plate_height = 0.3 + 0.22                    # [cm]
    hex_radius = float(hex_size*np.cos(30*(np.pi/180))) # [cm]
    fibre_separation = 0.51                      # [cm]
    
    vec_from_centre, host_vec, neighbour_vec, panel_name  = get_pmt_coordinates(host_hex, neighbour_hex, fname)

    # Points to the centre of the detctor!
    ######################################
    # vec_to_centre = -vec_from_centre
    ######################################
    # Rather than using the vector to the central PMT (comented out above)
    # try using the pointing vectors given in PANELINFO.ratdb - they agree
    # to within a few microns..... pretty sweet! 
    panel_table = ROOT.RAT.DB.Get().GetDefaultTable("PANELINFO", "")
    panel_numbers = panel_table.GetIArray("panel_number")
    panel_numbers_list = []
    for pan_num in panel_numbers:
        panel_numbers_list.append(pan_num)
    index =  panel_numbers_list.index(panel_name)
    vec_to_centre = ROOT.TVector3(panel_table.GetDArray("u")[index], panel_table.GetDArray("v")[index], panel_table.GetDArray("w")[index])
    #######################################
    vec_to_centre.SetMag(1.) 
    
    # Unit vector pointing from central hex to neighbour hex
    cent_to_neig_unit = neighbour_vec - vec_from_centre
    cent_to_neig_unit.SetMag(1.0)

    # Vector to mounted plate
    plate_vec = host_vec + cent_to_neig_unit*(hex_radius + plate_height)
    
    # Correct for 3D curve of PSUP
    curve_correction = (plate_vec - host_vec).Cross(vec_to_centre)
    curve_correction.SetMag(1.)

    # A & B positions
    positionA = plate_vec + curve_correction*fibre_separation
    positionB = plate_vec - curve_correction*fibre_separation

    # Some checks lifted stright from Sofia's code (docdb 1730)
    da=float(np.sin(positionA.Angle(vec_to_centre))*positionA.Mag())
    db=float(np.sin(positionB.Angle(vec_to_centre))*positionB.Mag())
    
    aa=cent_to_neig_unit.Angle(vec_to_centre)*(180./np.pi)
    ab=(positionA-positionB).Angle(vec_to_centre)*(180./np.pi)

    #if(da>500. or db>500. or np.abs(aa-90.)>0.5 or np.abs(ab-90.)>0.5):
    #    print "dist=", (positionA-positionB).Mag(), " in R=", positionA.Mag()-positionB.Mag()
    #    print "to center=", da, " and ",db
    #    print "Angle cor x dir=", aa, " and poA-poB x dir="

    return positionA, positionB

def get_pointing_angle(fibre):
    '''Get pointing angle of a fibre according to the info in the current
    database entry for that fibre
    '''
    try:
        entry = ROOT.RAT.DB.Get().GetLink("FIBRE", fibre)
        # Create vectors pointing to origin at centre of detector
        central_vector = ROOT.TVector3(entry.GetD("x"), entry.GetD("y"), entry.GetD("z"))
        pointing_vector = ROOT.TVector3(entry.GetD("u"), entry.GetD("v"), entry.GetD("w"))
    except Exception as e:
        print "Refernce to fibre %s does not exist in local ratdb" % fibre
        #raise e
    # Return the angle between the two vectors
    return 180 - central_vector.Angle(pointing_vector)*(180/np.pi)
        
def compare_position_calculations(fibres, pmt_hex, neighbour_hex, fname):
    '''Compare the x,y,z positions from our calculations with those in 
    the database
    '''
    count = 0
    for fibre in fibres:        

        posA, posB = calc_fibre_placement(pmt_hex[fibre], neighbour_hex[fibre], fname)
        table = ROOT.RAT.DB.Get().GetDefaultTable("FIBRE", fibre)
        try:
            table.GetD("x")
        except:
            print "###############################"
            print "Coundn't access x, y, z variables in database for fibre %s, skipping" % fibre
            continue        
        
        # Check for differences
        if fibre[-1] == "B":
            diff = table.GetD("x")/10. - posB.x()
        else:
            diff = table.GetD("x")/10. - posA.x()
        
        if np.abs(diff) > 1:
            print "###############################"
            print "%s:" % fibre
            print "Host: %s\nNeighbour: %s" % (pmt_hex[fibre], neighbour_hex[fibre])
            print "Database: \t%1.3f, %1.3f, %1.3f" % (table.GetD("x")/10., table.GetD("y")/10., table.GetD("z")/10.)
            if fibre[-1] == "B":
                print "Calculated :\t%1.3f, %1.3f, %1.3f" % (posB.x(), posB.y(), posB.z())
            else:
                print "Calculated :\t%1.3f, %1.3f, %1.3f" % (posA.x(), posA.y(), posA.z())

def make_new_db_files(fibres):
    '''Make new database files for jose
    '''
    for fibre in fibres:
        if fibre[:2] == "FT":
            angle = 0
        else:
            try:
                angle = int(round(get_pointing_angle(fibre), -1))
            except:
                "Fibre: %s does not exist in database" % fibre
                continue
        #print fibre, nodes[fibre], angle
        table = ROOT.RAT.DB.Get().GetDefaultTable("FIBRE", fibre)
        table.SetI("psup_pannel", nodes[fibre])
        table.SetI("pointing_angle", angle)
        table.SetI("fibre_status", 0)
        table.SaveAs("./new_tables/%s.json" % fibre)


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-i",dest="install",
                      help="Path to csv file detailing the current install status",
                      default="install_table 2016-05-09.csv")
    parser.add_option("-p",dest="pmt",
                      help="Path to csv file detailing the global co-ordinates of the PMTs",
                      default="PostRotation.txt")
    (options,args) = parser.parse_args()

    # Reset all root stuff
    ROOT.gROOT.Reset()

    # In in fibre install locations & hex cell positions
    nodes, fibres, pmt_hex, neighbour_hex = read_install_table(options.install)

    # Load rat defualts for fibre pos stuff
    ROOT.RAT.DB.Get().LoadDefaults()

    # make new database files for jose
    #make_new_db_files(fibres)

    # Compare position calculations
    compare_position_calculations(fibres, pmt_hex, neighbour_hex, options.pmt)
