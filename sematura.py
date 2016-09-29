#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

### Import required modules
import subprocess
import numpy as np
from sys import argv
from time import clock, time
import re
import glob
import sys
import os
import copy

### Import required modules from phenix software package
from dxtbx.format.Registry import Registry
from iotbx.phil import parse
from dxtbx.datablock import DataBlockFactory
from dials.array_family import flex
from dials.algorithms.indexing.real_space_grid_search \
     import indexer_real_space_grid_search as indexer

### Import user-defined parameters
from diffuse_params import *



class diffuse(object):

    def __init__(self, filenumber):
        if (os.path.isdir(work_dir+"/proc")) == True:
            pass
        else:
            subprocess.call(['mkdir', work_dir+"/proc"])
        if (os.path.isdir(work_dir+"/radial_averages")) == True:
            pass
        else:
            subprocess.call(['mkdir', work_dir+"/radial_averages"])
        if (os.path.isdir(work_dir+"/lattices")) == True:
            pass
        else:
            subprocess.call(['mkdir', work_dir+"/lattices"])

        self.name   = (work_dir+"/"+image_prefix+filenumber+".img")
        self.lunus  = (work_dir+"/proc/"+lunus_image_prefix
                       +filenumber+".img")
        self.radial = (work_dir+"/radial_averages/"+lunus_image_prefix
                       +filenumber+".asc")
        self.raw    = (raw_image_dir+"/"+image_prefix+filenumber+".cbf")
        self.latt   = (work_dir+"/lattices/"+diffuse_lattice_prefix
                       +"_diffuse_"+filenumber+".vtk")
        self.counts = (work_dir+"/lattices/"+diffuse_lattice_prefix
                       +"_counts_"+filenumber+".vtk")
        # self.ref = raw_image_dir+"/"+image_prefix+filenumber+".img"

    ### This function converts the raw data to the img format
    def cbf2img(self, filenumber):
        
        print "Converting file #"+filenumber+" to correct format"

        f = Registry.find(self.raw)
        img = f(self.raw)
        db = img.get_detectorbase()
        db.readHeader()
        db.read()
        db.show_header()
        print"Writing %s as %s"%(self.raw,self.name)

        db.debug_write(self.name)
        return

    ### This function uses Dials methods from Phenix to index the data
    def indexing(self):

        print target_cell,target_sg

        phil_scope_str='''
             include scope dials.algorithms.spot_finding.factory.phil_scope
             include scope dials.algorithms.indexing.indexer.index_only_phil_scope
             include scope dials.algorithms.refinement.refiner.phil_scope
             indexing.known_symmetry.unit_cell={0}
               .type = unit_cell
             indexing.known_symmetry.space_group={1}
               .type = space_group
             indexing.method=real_space_grid_search
             output.shoeboxes = False
            '''
        phil_scope = parse(phil_scope_str.format(target_cell,
                           target_sg), process_includes=True)
        params = phil_scope.extract()
        params.refinement.parameterisation.scan_varying = False
        params.spotfinder.filter.resolution_range = []
        params.indexing.scan_range = []
        params.refinement.parameterisation.crystal\
            .unit_cell.restraints.tie_to_target = []
        params.refinement.parameterisation.crystal\
            .unit_cell.restraints.tie_to_group = []

        filenames = [indexing_data_file_one, indexing_data_file_two,
                     indexing_data_file_three]


        datablock = DataBlockFactory.from_filenames(filenames)[0]

        observed = flex.reflection_table.from_observations(datablock, params)
        observed.as_pickle("strong.pickle")
        print "Number of observed reflections:", len(observed)

        working_params = copy.deepcopy(params)
        imagesets = datablock.extract_imagesets()


        print imagesets[0].get_beam()
        print imagesets[2].get_beam()
        print imagesets[0].get_beam() == imagesets[0].get_beam()
        print imagesets[1].get_beam() == imagesets[0].get_beam()
        print imagesets[2].get_beam() == imagesets[0].get_beam()

        print "indexing..."
        t0 = clock()

        idxr = indexer(observed, imagesets, params=working_params)

        tel = clock()-t0
        print "done indexing (",tel," sec)"

        indexed = idxr.refined_reflections
        experiments = idxr.refined_experiments
        print experiments.crystals()[0]

        global crystal_params

        crystal_params = experiments.crystals()[0]
        return
        
        ### don't keep...just to figure out how to extract variables
        # print type(experiments.crystals()[0])

        ### extract space group variable after indexing
        # print str(experiments.crystals()[0].get_space_group().info())

        ### extract unit cell variables after indexing
        # cella = experiments.crystals()[0].get_unit_cell().parameters()[0]
        # cellb = experiments.crystals()[0].get_unit_cell().parameters()[1]
        # cellc = experiments.crystals()[0].get_unit_cell().parameters()[2]
        # alpha = experiments.crystals()[0].get_unit_cell().parameters()[3]
        # beta  = experiments.crystals()[0].get_unit_cell().parameters()[4]
        # gamma = experiments.crystals()[0].get_unit_cell().parameters()[5]

    ### This function uses Lunus methods to remove the Bragg peaks from
    ### the raw data
    def debragg(self, filenumber):
        
        print "Removing beamstop shadow from image"
        subprocess.call(['punchim', self.name, punchim_xmin, punchim_xmax,
                          punchim_ymin, punchim_ymax, filenumber+"_tmp1.img"])
        print "Cleaning up edge of image"
        subprocess.call(['windim', filenumber+'_tmp1.img', windim_xmin,
                          windim_xmax, windim_ymin, windim_ymax,
                          filenumber+'_tmp2.img'])
        print "Setting detection threshold"
        subprocess.call(['thrshim', filenumber+'_tmp2.img', thrshim_min,
                          thrshim_max, filenumber+'_tmp3.img'])
        print "Correcting for beam polarization"
        subprocess.call(['polarim', filenumber+'_tmp3.img',
                         filenumber+'_tmp4.img', polarim_dist,
                         polarim_polarization, polarim_offset])
        print "Normalizing image"
        subprocess.call(['normim', filenumber+'_tmp4.img',
                         filenumber+'_tmp5.img', normim_tilt_x, normim_tilt_y])
        print "Removing Bragg peaks from image"
        subprocess.call(['modeim', filenumber+'_tmp5.img',
                         filenumber+'_tmp6.img', modeim_kernel_width,
                         modeim_bin_size])
        print "Cleaning up directory"
        subprocess.call(['cp', filenumber+'_tmp6.img', self.lunus])
        subprocess.call(['rm', filenumber+'_tmp1.img', filenumber+'_tmp2.img',
                         filenumber+'_tmp3.img',filenumber+'_tmp4.img',
                         filenumber+'_tmp5.img', filenumber+'_tmp6.img'])
        return


    ### This function uses Lunus methods to create a radially averaged copy
    ### of the image for scaling purposes
    def radial_avg(self, filenumber):

        print 'Averaging image'
        subprocess.call(['avgrim', self.lunus, filenumber+'_tmp.rf'])
        ### alternative is to keep the code below, make all input one string
        ### and set shell=True
        # subprocess.call(['binasc', '2', '<', 'tmp.rf', '>', self.radial])
        infile,outfile = filenumber+'_tmp.rf',self.radial
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = subprocess.Popen(
                    ['binasc', '2'],stdout=ouf,stdin=inf)
                proc.wait()
        return
    ### This function uses Lunus methods to create a radially averaged copy
    ### of the user-chosen reference image for scaling purposes
    def make_radial_ref(self):

        print 'Making reference statistic'
        ref_one = open(self.radial, 'r')
        ref_one_lines = ref_one.readlines()
        ref_one.close()
        # scale_inner_radius = scale_inner_radius - 1
        # this includes inner radius value
        ref_two = ref_one_lines[scale_inner_radius:scale_outer_radius]
        out_file = open('ref.asc', 'w')
        for line in ref_two:
            out_file.write(line)
        out_file.close()
        ### alternative is to keep the code below, make all input one string
        ### and set shell=True
        # subprocess.call(['binasc', '3', '<', 'ref.asc', '>', 'reference.rf'])
        infile,outfile = 'ref.asc', 'reference.rf'
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = subprocess.Popen(
                    ['binasc', '3'],stdout=ouf,stdin=inf)
                proc.wait()

        subprocess.call(['mulrf', 'reference.rf', 'reference.rf', 'xx.rf'])
        return

        # xx = subprocess.call(['avgrf', 'xx.rf'])

        # print "Cleaning up directory"

        # subprocess.call(['rm', 'xx.rf', 'ref.asc'])

        # return xx

    def scale_image(self, filenumber):
        print 'calculating scaling factor for image'

        ref_one = open(self.radial, 'r')
        ref_one_lines = ref_one.readlines()
        ref_one.close()
        # scale_inner_radius = scale_inner_radius - 1
        # this includes inner radius value
        ref_two = ref_one_lines[scale_inner_radius:scale_outer_radius]
        out_file = open(filenumber+'.asc', 'w')
        for line in ref_two:
            out_file.write(line)
        out_file.close()


        ### alternative is to keep the code below, make all input one string
        ### and set shell=True
        # subprocess.call(['binasc', '3', '<', filenum+'.asc',
        #                  '>', filenum+'.rf'])

        infile,outfile = filenumber+'.asc', filenumber+'.rf'
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = subprocess.Popen(
                    ['binasc', '3'],stdout=ouf,stdin=inf)
                proc.wait()

        subprocess.call(['mulrf', 'reference.rf', filenumber+'.rf',
                         filenumber+'_xy.rf'])
        subprocess.call(['mulrf', filenum+'.rf', filenumber+'.rf',
                         filenumber+'_yy.rf'])
        xx = float(subprocess.check_output(['avgrf', 'xx.rf']))
        xy = float(subprocess.check_output(['avgrf', filenumber+'_xy.rf']))
        yy = float(subprocess.check_output(['avgrf', filenumber+'_yy.rf']))

        global scale_factor, scale_factor_error
        scale_factor = xx / xy
        scale_factor_error = np.sqrt(xx+yy*scale_factor*scale_factor
                                     -2.*scale_factor*xy)/np.sqrt(xx)

        # return scale_factor, scale_factor_error

        print "This image has a scale factor of "+str(scale_factor)
        print "This image has a scale factor error of "+str(scale_factor_error)
        return

    # this function maps diffuse data to a 3D lattice
    def procimg(self, Isize1,Isize2,scale_factor,mask_tag,A_matrix,rvec,
                DATA,latxdim,latydim,latzdim):
        
        global lat

        tmid = clock()

        # define the lattice indices at which h,k,l = 0,0,0
        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1
        # total number of voxels in the lattice
        latsize = latxdim*latydim*latzdim
        # generate lattice to store data & counts for averaging
        lat = np.zeros(latsize*2, dtype=np.float32).reshape((2,latsize))

        # Fetch A_matrix from processed image using Dials
        alli = np.asanyarray(A_matrix)
        alli = np.reshape(alli, (3,3))
        # calculate h,k,l for every data point
        H = np.tensordot(alli,rvec,(1,1))
        H = np.transpose(H)

        tel = str(clock()-tmid)
        print "Calculated h,k,l for ",len(H[:]),\
              " data points ("+tel+" seconds)"
        

        # fetch data from processed image using Dials
        val = np.asanyarray(DATA, dtype=int)

        # rearrange order of data points to match h,k,l matrix (H)
        val.shape = (Isize2,Isize1)
        val = np.transpose(val)
        val = val.flatten()

        # isolate all h's, k's, and l's
        ii = H[:,0]
        jj = H[:,1]
        kk = H[:,2]

        # adjust hkls
        H[H>=0]+=0.5
        H[H<0]-=0.5
        H = np.asanyarray(H, dtype=int)

        # adjusted hkls
        i = H[:,0]
        j = H[:,1]
        k = H[:,2]

        ### calculate the displacement of this data point
        ### from the nearest Miller index
        # make a mask that eliminates diffuse points
        # that are too close to bravais lattice (h)
        dimask = abs(ii-i)
        dimask[dimask<0.25]=0
        dimask[dimask!=0]=1.0
        # make a mask that eliminates diffuse points 
        # that are too close to bravais lattice (k)
        djmask = abs(jj-j)
        djmask[djmask<0.25]=0
        djmask[djmask!=0]=1.0
        # make a mask that eliminates diffuse points
        # that are too close to bravais lattice (l)
        dkmask = abs(kk-k)
        dkmask[dkmask<0.25]=0
        dkmask[dkmask!=0]=1.0

        # adjust hkls to map onto diffuse lattice
        i = i + i0
        j = j + j0
        k = k + k0

        # make an array of indices to map data onto the diffuse lattice
        index = k*latxdim*latydim + j*latxdim + i
        index = np.asanyarray(index, dtype=int)

        # create a mask to eliminate any data point outside 1 unit cell (h)
        imask = i
        imask[imask<0]=0
        imask[imask>=latxdim]=0
        imask[imask!=0]=1.0
        # create a mask to eliminate any data point outside 1 unit cell (k)
        jmask = j
        jmask[jmask<0]=0
        jmask[jmask>=latydim]=0
        jmask[jmask!=0]=1.0
        # create a mask to eliminate any data point outside 1 unit cell (l)
        kmask = k
        kmask[kmask<0]=0
        kmask[kmask>=latzdim]=0
        kmask[kmask!=0]=1.0

        # eliminate data with negative intensities
        val[val<0]=0.0
        # eliminate data marked with mask_tag by lunus software
        val[val>=mask_tag]=0.0
        # apply masks generated above
        val = val*imask*jmask*kmask*dimask*djmask*dkmask

        # map the data onto the diffuse lattice using the indices created above
        np.add.at(lat[0], index, val)
        # keep track of the number of data points
        # added at each lattice point (for averaging)
        val[val!=0]=1
        np.add.at(lat[1], index, val)

        return lat

    # This function integrates the diffuse scattering using the procimg fxn
    def integrator(self):

        global lt, ct, latxdim, latydim, latzdim, i0, j0, k0, latsize

        res = float(resolution)
        latxdim = (int(cella/res)+1)*2
        latydim = (int(cellb/res)+1)*2
        latzdim = (int(cellc/res)+1)*2


        latsize = latxdim*latydim*latzdim
        print "Lattice size = ",latsize
        lt = np.zeros(latsize, dtype=np.float32)
        ct = np.zeros(latsize, dtype=np.float32)



        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1
        mask_tag = 32767



        imgname = self.lunus

        print "processing file %s with scale factor %f"%(imgname,scale_factor)

        import dxtbx
        img = dxtbx.load(imgname)
        detector = img.get_detector()
        print detector
        beam = img.get_beam()
        scan = img.get_scan()
        gonio = img.get_goniometer()

        print "transform pixel numbers to mm positions and rotational degrees"


        print "Creating pixel map..."
        t0 = clock()


        lab_coordinates = flex.vec3_double()
        for panel in detector: 
            pixels = flex.vec2_double(panel.get_image_size())
            mms = panel.pixel_to_millimeter(pixels)
            lab_coordinates.extend(panel.get_lab_coord(mms))

        # generate s1 vectors
        s1_vectors = lab_coordinates.each_normalize() * (1/beam.get_wavelength())
        # Generate x vectors
        x_vectors = s1_vectors - beam.get_s0()

        print "there are ",x_vectors.size()," elements in x_vectors"
        tel = clock()-t0
        print "done creating pixel map (",tel," sec)"

        print "transform to laboratory axis reciprocal space coordinates"

        print "transform to fractional miller indices and populate diffuse lattice"

        crystal = copy.deepcopy(crystal_params)
        axis = gonio.get_rotation_axis()
        start_angle, delta_angle = scan.get_oscillation()
        crystal.rotate_around_origin(axis, start_angle + (delta_angle/2), deg=True)
        A_matrix = crystal.get_A().inverse()


        telmatmul=0
        t0 = clock()
        latit = None
        for panel_id, panel in enumerate(detector):

            Isize1, Isize2 = panel.get_image_size()
            print "Isize1 = ",Isize1,", Isize2 = ",Isize2
            print "there are ",Isize1*Isize2," pixels in this diffraction image"
            if len(detector) > 1:
                DATA = img.get_raw_data(panel_id)
            else:
                DATA = img.get_raw_data()

            tmp_latit = self.procimg(Isize1,Isize2,scale_factor,mask_tag,A_matrix,
                                x_vectors,DATA,latxdim,latydim,latzdim)
            if latit is None:
                latit = tmp_latit
            else:
                latit += tmp_latit
        tel = clock()-t0
        print "done integrating diffuse scattering (",tel," sec wall clock time)"
        t0 = clock()
        # accumulate integration data into a single lattice

        lt = np.add(lt,latit[0])
        ct = np.add(ct,latit[1])
        tel = clock()-t0
        print "Took ",tel," secs to update the lattice"

        return lt, ct

    # This function writes a vtk file containing the diffuse lattice
    def latout(self):

        make1 = self.latt
        # write lattice to output file
        vtkfile = open(make1,"w")

        a_recip = 1./cella
        b_recip = 1./cellb
        c_recip = 1./cellc

        print >>vtkfile,"# vtk DataFile Version 2.0"
        print >>vtkfile,"lattice_type=PR;unit_cell={0};space_group={1};".format(target_cell,target_sg)
        print >>vtkfile,"ASCII"
        print >>vtkfile,"DATASET STRUCTURED_POINTS"
        print >>vtkfile,"DIMENSIONS %d %d %d"%(latxdim,latydim,latzdim)
        print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
        print >>vtkfile,"ORIGIN %f %f %f"%(-i0*a_recip,-j0*b_recip,-k0*c_recip)
        print >>vtkfile,"POINT_DATA %d"%(latsize)
        print >>vtkfile,"SCALARS volume_scalars float 1"
        print >>vtkfile,"LOOKUP_TABLE default\n"

        index = 0
        for k in range(0,latzdim):
            for j in range(0,latydim):
                for i in range(0,latxdim):
                    print >>vtkfile,lt[index],
                    index += 1
                print >>vtkfile,""

        vtkfile.close()
        return
    
    # This function writes a vtk file containing the count
    # of data points contributing to each point on the diffuse lattice
    def ctout(self):
        make2 = self.counts

        vtkfile = open(make2,"w")

        a_recip = 1./cella
        b_recip = 1./cellb
        c_recip = 1./cellc

        print >>vtkfile,"# vtk DataFile Version 2.0"
        print >>vtkfile,"lattice_type=PR;unit_cell={0};space_group={1};".format(target_cell,target_sg)
        print >>vtkfile,"ASCII"
        print >>vtkfile,"DATASET STRUCTURED_POINTS"
        print >>vtkfile,"DIMENSIONS %d %d %d"%(latxdim,latydim,latzdim)
        print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
        print >>vtkfile,"ORIGIN %f %f %f"%(-i0*a_recip,-j0*b_recip,-k0*c_recip)
        print >>vtkfile,"POINT_DATA %d"%(latsize)
        print >>vtkfile,"SCALARS volume_scalars float 1"
        print >>vtkfile,"LOOKUP_TABLE default\n"

        index = 0
        for k in range(0,latzdim):
            for j in range(0,latydim):
                for i in range(0,latxdim):
                    print >>vtkfile,ct[index],
                    index += 1
                print >>vtkfile,""

        vtkfile.close()
        return

    # This function calculates a mean lattice using all diffuse and counts
    # files
    def mean_lattice(self):
        
        # i = 1
        # num_imgs = 360 ###!!!put this in params or define using iglob!!!###

        # # initialize lattice before summing
        # latnumber = '{:05.0f}'.format(i)
        lattice1   = (work_dir+"/lattices/"+diffuse_lattice_prefix
                      +"_diffuse_"+filenum+".vtk")
        

        subprocess.call(['vtk2lat', lattice1, 'temp.lat'])
        subprocess.call(['constlt', 'temp.lat', 'temp_diffuse_sum.lat', '0.0'])
        subprocess.call(['constlt', 'temp.lat', 'temp_counts_sum.lat', '0.0'])

        # sum diffuse lattices and counts
        for item in diffuse_file_list:
            # latnumber = '{:05.0f}'.format(i)
            # lattice   = (raw_image_dir+"/lattices/"+diffuse_lattice_prefix
            #              +"_diffuse_"+latnumber+".vtk")
            # counts   = (raw_image_dir+"/lattices/"+diffuse_lattice_prefix
                         # +"_counts_"+latnumber+".vtk")

            subprocess.call(['vtk2lat', item, 'temp_diffuse.lat'])
            
            subprocess.call(['sumlt', 'temp_diffuse_sum.lat',
                             'temp_diffuse.lat', 'temp.lat'])
            subprocess.call(['mv', 'temp.lat', 'temp_diffuse_sum.lat'])
            
            # i+=1
        for item in counts_file_list:
            subprocess.call(['vtk2lat', item, 'temp_counts.lat'])
            subprocess.call(['sumlt', 'temp_counts_sum.lat',
                             'temp_counts.lat', 'temp.lat'])
            subprocess.call(['mv', 'temp.lat', 'temp_counts_sum.lat'])

        
        subprocess.call(['divlt', 'temp_diffuse_sum.lat', 'temp_counts_sum.lat'
                         , 'temp_mean.lat'])
        subprocess.call(['lat2vtk', 'temp_mean.lat', mean_lattice_file])
        subprocess.call(['rm', 'temp_diffuse_sum.lat', 'temp_counts_sum.lat'])
        subprocess.call(['mv', 'temp_mean.lat', mean_lat_file])
        return


    def symmetrize(self):

        sg_conv = open(lunus_dir+"/analysis/sg_pg_lunus.csv","r")
        sg_conv_lines = sg_conv.readlines()

        lines = None
        lunus_key = ""
        laue_class = ""
        for line in sg_conv_lines:
            line = line.split('\r')
            lines = line


        for item in lines:
            item = item.split(',')
            if item[0] == target_sg:
                lunus_key += item[4]
                laue_class += item[3]
        
        print lunus_key
        print laue_class
        # sg_conv = pd.read_csv(work_dir+"/analysis/sg_pg_lunus.csv", low_memory=False)
        # lunus_key = sg_conv['lunus_number'][target_sg]

        subprocess.call(['xflt', mean_lat_file, 'tmp_xf.lat', '1'])
        print "Symmetrizing lattice based on symmetry operators for Laue Class: "+laue_class
        subprocess.call(['symlt', 'tmp_xf.lat', sym_lat_file, lunus_key])
        subprocess.call(['lat2vtk', sym_lat_file, sym_lattice_file])
        return


###----------------Begin Variable Definitions-------------------###

lunus_dir = subprocess.check_output(['which', 'symlt'])
lunus_dir = lunus_dir.replace("/c/bin/symlt\n","")
work_dir = subprocess.check_output("pwd")
work_dir = work_dir.replace("\n","")

### variables to be plugged in via command line
script, filenum = argv
filenum = float(filenum)
filenum = '{:05.0f}'.format(filenum)

### match number formatting for reference image
reference_image_number = float(reference_image_number)
reference_image_number = '{:05.0f}'.format(reference_image_number)

raw_file_list = glob.glob(raw_image_dir+"/"+image_prefix+"*.cbf")
lunus_files = glob.iglob(work_dir+"/proc/"+lunus_image_prefix+"*.img")
radial_files = glob.iglob(work_dir+"/radial_averages/"+lunus_image_prefix
                          +"*.asc")





###-----------------End Variable Definitions--------------------###


###---------Runs fxns within Diffuse class...all work------------###

data = diffuse(filenum)
reference = diffuse(reference_image_number)
data.cbf2img(filenum)
# make lunus processed image for file of interest
data.debragg(filenum)

# # make reference image for scaling
### generate list of files (don't use)

# files = subprocess.check_call['ls', './*.sh']
# files = subprocess.check_output('ls /Users/student/Desktop/vanben_diffuse_analysis/*.sh', shell=True)
# print files

### generate list of files (use)



# for item in files:
#     print item # ['file1.bc', 'file2.bc']


if reference.lunus in lunus_files:
    pass
else:
    reference.cbf2img(reference_image_number)
    reference.debragg(reference_image_number)




data.radial_avg(filenum)

# make reference image for scaling
if reference.radial in radial_files:
    pass
else:
    reference.radial_avg(reference_image_number)

reference.make_radial_ref()
data.scale_image(filenum)
data.indexing()
data.integrator()
data.latout()
data.ctout()

lattice_file_list = glob.glob(work_dir+"/lattices/"+diffuse_lattice_prefix
                              +"*.vtk")
diffuse_file_list = glob.glob(work_dir+"/lattices/"+diffuse_lattice_prefix
                              +"_diffuse_*.vtk")
counts_file_list = glob.glob(work_dir+"/lattices/"+diffuse_lattice_prefix
                              +"_counts_*.vtk")

mean_lat_file   = (work_dir+"/lattices/"+diffuse_lattice_prefix
                       +"_mean.lat")
mean_lattice_file = (work_dir+"/lattices/"+diffuse_lattice_prefix
                       +"_mean.vtk")
sym_lat_file = (work_dir+"/lattices/"+diffuse_lattice_prefix
                      +"_mean_sym.lat")
sym_lattice_file = (work_dir+"/lattices/"+diffuse_lattice_prefix
                      +"_mean_sym.vtk")

# calculate mean lattice only once all images are processed
if (len(diffuse_file_list)+len(counts_file_list)) < (2*len(raw_file_list)):
    pass
elif (len(diffuse_file_list)+len(counts_file_list)) == (2*len(raw_file_list)):
    data.mean_lattice()
else:
    print "Mean lattice already exists"



# calculate symmetrized lattice only if mean lattice exists
if mean_lattice_file in lattice_file_list:
    # data.sg_pg_lunus()
    data.symmetrize()
else:
    pass
###---------End of run within Diffuse class...all work------------###











