#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

"""
The purpose of this library is to develop automated methods
for extracting the diffuse X-ray scattering from raw
diffraction frames.

Author(s):
Alexander M. Wolff      (primary author)
Michael E. Wall         (primary author)
Andrew Van Benschoten   (contributed initial methods)

"""

### Import required functions
import numpy as np
from subprocess import call, check_output, Popen
from glob import glob
from sys import argv
from copy import deepcopy
import time

from libtbx.phil import parse
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from os.path import basename, isdir
from os import getcwd
from dxtbx import load
from scitbx.matrix import sqr, col
from dials.array_family import flex
from lunus import LunusDIFFIMAGE, LunusLAT3D
from dxtbx.format.FormatCBFMini import FormatCBFMini




lunus_sym_lib = {'-1': 0, '2/m': -1, 'mmm': -2, '4/m': -3, '4/mmm': -4, \
'-3': -5, '-3m': -6, '6/m': -7, '6/mmm': -8, 'm-3': -9, 'm-3m': -10}

# lunus_dir = check_output(['which', 'symlt'])
# lunus_dir = lunus_dir.replace("/c/bin/symlt\n","")
work_dir = getcwd()
# work_dir = check_output("pwd")
# work_dir = work_dir.replace("\n","")

class DiffuseExperiment(object):

    def __init__(self):

        if (isdir(work_dir+"/lattices")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/lattices"])
        if (isdir(work_dir+"/arrays")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/arrays"])


    def read_experiment_file(self, experiment_file):

        ### open DIALS json file
        phil_scope_str='''
            experiments = 'example_refined_experiments.json'
              '''
        phil_scope = parse(phil_scope_str, process_includes=True)
        parser = OptionParser(
            phil=phil_scope,
            check_format=False,
            read_experiments=True)
        params, options = parser.parse_args(args=[experiment_file], show_diff_phil=True)
        experiments = flatten_experiments(params.input.experiments)
        exp_xtal = experiments.crystals()[0]

        ### define useful attributes
        self.crystal = exp_xtal
        uc = self.crystal.get_unit_cell()
        uc_nospace = str(uc).replace(" ", "")
        uc_nospace_noparen = uc_nospace[1:-1]
        self.unit_cell = uc_nospace_noparen
        self.space_group = self.crystal.get_space_group()
        self.laue_group = self.space_group.laue_group_type()
        # self.a_matrix = crystal.get_A()
        self.experiments = experiments


    def read_reflection_file(self, reflection_file):

        ### open DIALS pickle file
        phil_scope_str='''
            reflections = 'example_refined.pickle'
              '''
        phil_scope = parse(phil_scope_str, process_includes=True)
        parser = OptionParser(
            phil=phil_scope,
            check_format=False,
            read_reflections=True)
        params, options = parser.parse_args(args=[reflection_file], show_diff_phil=True)
        reflections = flatten_reflections(params.input.reflections)
        self.reflections = reflections


    def xds_to_dials(exp_dir):
        call(["dials.import_xds", "xds_file=/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/integrate/XPARM.XDS", "/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/integrate/", "filename=/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/refine/experiments.json"])
        call(["dials.import_xds", "/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/integrate/INTEGRATE.HKL", "/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/refine/experiments.json", "method=reflections", "filename=/Volumes/beryllium/sbgrid_project/bragg_processing/110/3dii/DEFAULT/NATIVE/SWEEP1/refine/second.pickle"])

        return



class DiffuseImage(DiffuseExperiment):

    def __init__(self, filename):

        self.raw    = filename
        try:
            left, right = self.raw.split(".")
            self.lunus = left+"_lunus."+right
        except:
            left, middle, right = self.raw.split(".")
            self.lunus = left+"_lunus."+middle+"."+right
        base = basename(self.raw)
        try:
            self.id, self.filetype = base.split(".")
        except:
            self.id, self.filetype, self.compression = base.split(".")

        self.array = (work_dir+"/arrays/"+self.id+".npz")

        print("raw file is: {}".format(filename))
        print("lunus processed file is: {}".format(self.lunus))
        print(self.id)
        print(self.filetype)
        print(self.array)


    def set_general_variables(self):

        self.frame = load(self.raw)
        self.detector = self.frame.get_detector()
        self.beam = self.frame.get_beam()
        self.s0 = self.beam.get_s0()
        self.gonio = self.frame.get_goniometer()
        self.scan = self.frame.get_scan()

        self.lab_coordinates = flex.vec3_double()
        for panel in self.detector:
            self.beam_center_mm_x, self.beam_center_mm_y = col(panel.get_beam_centre(self.s0))
            pixels = flex.vec2_double(panel.get_image_size())
            mms = panel.pixel_to_millimeter(pixels)
            self.lab_coordinates.extend(panel.get_lab_coord(mms))
            self.Isizex, self.Isizey = panel.get_image_size()
            self.beam_center_x, self.beam_center_y = col(panel.get_beam_centre_px(self.s0))
            self.detector_distance = panel.get_distance()
            thrshim_min, thrshim_max = panel.get_trusted_range()
            self.pixel_size = panel.get_pixel_size()[0]

        self.raw_data = self.frame.get_raw_data()

        if thrshim_min < 0 :
            self.thrshim_min = int(0)
        else:
            self.thrshim_min = thrshim_min

        if thrshim_max > 32767:
            self.thrshim_max = int(32767)
        else:
            self.thrshim_max = int(thrshim_max)

        self.polarization_fraction = self.beam.get_polarization_fraction()
        self.polarization_offset = 0.0
        self.cassette_x = 0.0
        self.cassette_y = 0.0
        self.windim_xmax     = int(self.Isizex)-100          # right border for processed image (pixels)
        self.windim_xmin     = 100          # left border for processed image (pixels)
        self.windim_ymax     = int(self.Isizey)-100     # top border for processed image (pixels)
        self.windim_ymin     = 100          # bottom border for processed image (pixels)
        ### beamstop borders
        self.punchim_xmax    = int(self.Isizex)          # right border of beam stop shadow (pixels)
        self.punchim_xmin    = int(self.beam_center_x)-80          # left border of beam stop shadow (pixels)
        self.punchim_ymax    = int(self.beam_center_y)+100          # top border of beam stop shadow (pixels)
        self.punchim_ymin    = int(self.beam_center_y)-40          # bottom border of beam stop shadow (pixels)
        self.mode_filter_footprint = int(20)

        return


    def remove_bragg_peaks(self, reference=None, write_lunus=None):

        start_cpu = time.clock()
        start_real = time.time()
        image = LunusDIFFIMAGE()
        image.set_image(self.raw_data)
        print 'Removing beamstop shadow from image'
        image.LunusPunchim(self.punchim_xmin, self.punchim_ymin, self.punchim_xmax, self.punchim_ymax)
        print 'Cleaning up edge of image'
        image.LunusWindim(self.windim_xmin, self.windim_ymin, self.windim_xmax, self.windim_ymax)
        print "Setting detection threshold"
        image.LunusThrshim(self.thrshim_min, self.thrshim_max)
        print "Removing Bragg peaks from image"
        image.LunusModeim(self.mode_filter_footprint)
        print "Mode filter finished."
        print "Correcting for beam polarization"
        image.LunusPolarim(self.beam_center_mm_x, self.beam_center_mm_y, self.detector_distance, self.polarization_fraction, self.polarization_offset, self.pixel_size)
        print "Solid-angle normalization of image"
        image.LunusNormim(self.beam_center_mm_x, self.beam_center_mm_y, self.detector_distance, self.cassette_x, self.cassette_y, self.pixel_size)
        
        self.lunus_data_scitbx = image.get_image();

        self.radial_avg = image.LunusRadialAvgim(self.beam_center_mm_x, self.beam_center_mm_y, self.pixel_size)
        self.lunus_data = self.lunus_data_scitbx.as_numpy_array()
        print("shape of lunus data at step 1: {}".format(self.lunus_data.shape))


        end_cpu = time.clock()
        end_real = time.time()
        print("%f Real Seconds" % (end_real - start_real))
        print("%f CPU seconds" % (end_cpu - start_cpu))


        ### write files, if desired
        if write_lunus:
            FormatCBFMini.as_file(self.detector,self.beam,self.gonio,self.scan,self.lunus_data_scitbx,self.lunus)
        else:
            pass

        if reference:
            np.savez("reference_radial_average", rad=self.radial_avg)
        else:
            pass

        return


    def radial_average(self, reference=None):
        image = LunusDIFFIMAGE()
        working_img = self.lunus_data_scitbx.reshape(flex.grid(self.Isizey, self.Isizex))
        image.set_image(working_img)
        self.radial_avg = image.LunusRadialAvgim(self.beam_center_mm_x, self.beam_center_mm_y, self.pixel_size)

        if reference:
            np.savez("reference_radial_average", rad=self.radial_avg)
        else:
            pass
        
        return
    # def radial_average(self, reference=None):

    #     y, x = np.indices((self.lunus_data.shape))
    #     r = np.sqrt((x - self.beam_center_x)**2 + (y - self.beam_center_y)**2)
    #     r = r.astype(np.int)
    #     data_mask = np.array(self.lunus_data, dtype=bool)
    #     data_mask[self.lunus_data<1e-6]=False
    #     r=r[data_mask]
    #     data=self.lunus_data[data_mask]
    #     tbin = np.bincount(r.ravel(), data.ravel())
    #     nr = np.bincount(r.ravel())
    #     with np.errstate(invalid='ignore', divide='ignore'):
    #         self.radial_avg = tbin / nr

    #     if reference:
    #         np.savez("reference_radial_average", rad=self.radial_avg)
    #     else:
    #         pass

    #     return


    def scale_factor(self):

        ref = np.load("reference_radial_average.npz")
        ref_rad = ref["rad"]
        xx = np.average(np.multiply(ref_rad[100:800],ref_rad[100:800]))
        xy = np.average(np.multiply(ref_rad[100:800],self.radial_avg[100:800]))
        self.scale_factor = xx / xy
        print "This image has a scale factor of "+str(self.scale_factor)

        return


    def crystal_geometry(self, xtal):

        self.crystal = xtal
        uc = self.crystal.get_unit_cell()
        uc_nospace = str(uc).replace(" ", "")
        uc_nospace_noparen = uc_nospace[1:-1]
        self.unit_cell = uc_nospace_noparen
        self.space_group = self.crystal.get_space_group()
        self.laue_group = self.space_group.laue_group_type()
        self.cella, self.cellb, self.cellc, self.alpha, self.beta, self.gamma = uc.parameters()

        s1_vectors = self.lab_coordinates.each_normalize() * (1/self.beam.get_wavelength())
        self.r_vectors = s1_vectors - self.beam.get_s0()

        crystal = deepcopy(xtal)
        axis = self.gonio.get_rotation_axis()
        start_angle, delta_angle = self.scan.get_oscillation()
        crystal.rotate_around_origin(axis, start_angle + (delta_angle/2), deg=True)
        self.A_matrix = crystal.get_A()

        return


    def setup_diffuse_lattice(self, d_min):

        ### set resolution of diffuse lattice
        self.latxdim = (int(self.cella/d_min)+1)*2
        self.latydim = (int(self.cellb/d_min)+1)*2
        self.latzdim = (int(self.cellc/d_min)+1)*2
        self.latsize = self.latxdim*self.latydim*self.latzdim
        print("Lattice size = {}".format(self.latsize))
        self.lt = np.zeros(self.latsize, dtype=np.float32)
        self.ct = np.zeros(self.latsize, dtype=np.float32)
        ### set origin and tag value for null data
        self.i_0=self.latxdim/2-1
        self.j_0=self.latydim/2-1
        self.k_0=self.latzdim/2-1
        self.mask_tag = 32767

        return


    def image_to_lattice(self):
        
        # global lat
        tmid = time.clock()

        ### generate lattice to store data & counts for averaging
        lat = np.zeros(self.latsize*2, dtype=np.float32).reshape((2,self.latzdim,self.latydim,self.latxdim))

        ### fetch A_matrix from Dials output after indexing
        a_mat = np.array(self.A_matrix)
        a_mat.shape = (3,3)
        Amat = np.asmatrix(a_mat)
        from numpy.linalg import inv
        a_inv = inv(Amat)
        # a_inv = a = np.array([[42.7267962919, -4.43572592458, -2.81263484758],[-6.13535435886, -50.434871775, -13.6629067462],[-3.20694000173, 23.7226323386, -86.1289788991]])

        ### calculate h,k,l for every data point
        print "Calculating h,k,l for each data point from image"
        H = np.tensordot(a_inv, self.r_vectors, (1,1))
        H = np.transpose(H)
        tel = str(time.clock()-tmid)
        print "Calculated h,k,l for ",len(H[:]),\
              " data points ("+tel+" seconds)"

        ### fetch data from Lunus processed image using Dials methods
        val = np.array(self.lunus_data_scitbx, dtype=int)
        print("shape of lunus data at step 2: {}".format(val.shape))


        ### rearrange order of data points to match h,k,l matrix (H)
        val.shape = (self.Isizey,self.Isizex)
        print("shape of lunus data at step 3: {}".format(val.shape))

        val = np.transpose(val)
        print("shape of lunus data at step 4: {}".format(val.shape))

        val = val.flatten()
        print("shape of lunus data at step 5: {}".format(val.shape))


        ### isolate all h's, k's, and l's
        i = H[:,0]
        j = H[:,1]
        k = H[:,2]

        ### adjust hkls to integer values
        H_int = np.copy(H)
        H_int[H_int>=0]+=0.5
        H_int[H_int<0]-=0.5
        H_int = np.array(H_int, dtype=int)

        ### isolate all integer h's, k's, and l's
        i_int = H_int[:,0]
        j_int = H_int[:,1]
        k_int = H_int[:,2]

        ### calculate the displacement of every data point
        ### from the nearest Miller index
        delta_i = abs(i-i_int)
        delta_j = abs(j-j_int)
        delta_k = abs(k-k_int)

        ### adjust coordinates so that origin is centered on diffuse lattice
        i_int += int(self.i_0)
        j_int += int(self.j_0)
        k_int += int(self.k_0)

        ### make a boolean mask to be used for data & hkls
        data_mask = np.array(val, dtype=bool)
        ### mask data that are too close to bravais lattice
        data_mask[delta_i<0.25] = False
        data_mask[delta_j<0.25] = False
        data_mask[delta_k<0.25] = False
        ### mask data outside unit cell on diffuse lattice
        data_mask[i_int<0] = False
        data_mask[i_int>=self.latxdim] = False
        data_mask[j_int<0] = False
        data_mask[j_int>=self.latydim] = False
        data_mask[k_int<0] = False
        data_mask[k_int>=self.latzdim] = False
        ### mask data with negative intensities
        data_mask[val<=0] = False
        ### mask data marked with mask_tag by lunus software
        data_mask[val>=self.mask_tag] = False
        
        ### apply mask to data and associated hkls
        val = val[data_mask]
        i_int = i_int[data_mask]
        j_int = j_int[data_mask]
        k_int = k_int[data_mask]

        ### map the data onto the diffuse lattice
        np.add.at(lat[0], [k_int,j_int,i_int], val)
        diff = lat[0].flatten()
        diff_scaled = np.multiply(diff,self.scale_factor)
        ### keep track of the number of data points
        ### added at each lattice point (for averaging)
        val[val!=0]=1
        np.add.at(lat[1], [k_int,j_int,i_int], val)
        coun = lat[1].flatten()

        # diff_scaled[coun==0] = 32767

        lat = [diff_scaled,coun]

        return lat


    def integrate_diffuse(self):

        t0 = time.clock()
        latit = None
        for panel_id, panel in enumerate(self.detector):

            print "there are ",self.Isizex * self.Isizey," pixels in this diffraction image"
            # if len(detector) > 1:
            #     DATA = img.get_raw_data(panel_id)
            # else:
            #     DATA = img.get_raw_data()
            tmp_latit = self.image_to_lattice()
            if latit is None:
                latit = tmp_latit
            else:
                latit += tmp_latit

        tel = time.clock()-t0
        print("Integrating diffuse scattering took {} sec".format(tel))
        lt = np.add(self.lt,latit[0])
        ct = np.add(self.ct,latit[1])
        np.savez(self.array, lt=lt, ct=ct)

        return


    def array_to_vtk(self):

        in_file = (work_dir+"/arrays/"+self.id+"_mean.npz")
        key = "mean_lt"
        array = np.load(in_file)
        array_values = array[key]
        out_file = in_file.replace("/arrays/","/lattices/")
        out_file = out_file.replace(".npz",".vtk")
        vtkfile = open(out_file,"w")

        array_values = array_values.flatten()

        a_recip = 1./self.cella
        b_recip = 1./self.cellb
        c_recip = 1./self.cellc

        print >>vtkfile,"# vtk DataFile Version 2.0"
        print >>vtkfile,"lattice_type=PR;unit_cell={0},{1},{2},{3},{4},{5};space_group={6};".format(self.cella,self.cellb,self.cellc,self.alpha,self.beta,self.gamma, self.unit_cell)
        print >>vtkfile,"ASCII"
        print >>vtkfile,"DATASET STRUCTURED_POINTS"
        print >>vtkfile,"DIMENSIONS %d %d %d"%(self.latxdim,self.latydim,self.latzdim)
        print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
        print >>vtkfile,"ORIGIN %f %f %f"%(-self.i_0*a_recip,-self.j_0*b_recip,-self.k_0*c_recip)
        print >>vtkfile,"POINT_DATA %d"%(self.latsize)
        print >>vtkfile,"SCALARS volume_scalars float 1"
        print >>vtkfile,"LOOKUP_TABLE default\n"

        index = 0
        for k in range(0,self.latzdim):
            for j in range(0,self.latydim):
                for i in range(0,self.latxdim):
                    print >>vtkfile,array_values[index],
                    index += 1
                print >>vtkfile,""

        vtkfile.close()
        self.vtk = out_file
        return


    def mean_lattice(self):

        print "Calculating mean lattice"
        intensity_file_list = glob(work_dir+"/arrays/*.npz")
        sum_file = np.load(intensity_file_list[0])
        sum_int = sum_file['lt']
        sum_ct = sum_file['ct']
        sum_int[:]=0
        sum_ct[:]=0
        for item in intensity_file_list:
            data = np.load(item)
            sum_int += data['lt']
            sum_ct += data['ct']
            data.close()
        sum_int_masked = np.ma.array(sum_int,mask=sum_ct==0)
        sum_ct_masked = np.ma.array(sum_ct,mask=sum_ct==0)
        mean_lt = np.ma.divide(sum_int_masked,sum_ct_masked)
        mean_lt.set_fill_value(-32768)
        self.mean_lattice = mean_lt.filled()
        out = (work_dir+"/arrays/"+self.id+"_mean.npz")
        np.savez(out, mean_lt=self.mean_lattice)

        return


    def average_symmetry_mates(self, lat=None, vtk=None):

        ### filenames for outputs
        mean_vtk_file   = (work_dir+"/lattices/"+self.id+"_mean.vtk")
        mean_lat_file   = (work_dir+"/lattices/"+self.id+"_mean.lat")
        sym_lat_file    = (work_dir+"/lattices/"+self.id+"_mean_sym.lat")
        sym_vtk_file    = (work_dir+"/lattices/"+self.id+"_mean_sym.vtk")
        aniso_lat_file  = (work_dir+"/lattices/"+self.id+"_mean_sym_aniso.lat")
        aniso_vtk_file  = (work_dir+"/lattices/"+self.id+"_mean_sym_aniso.vtk")

        lat = LunusLAT3D()
        ### numpy arrays read into flex arrays don't store data properly
        # working_lattice = flex.int([np.int(i) for i in self.mean_lattice.flatten()])
        # working_lattice.reshape(flex.grid(self.latzdim,self.latydim,self.latxdim))
        lat.LunusReadvtk(self.vtk)
        # lat.set_lattice(working_lattice)
        if lat:
            lat.LunusWritelt(mean_lat_file)
        else:
            pass 
        if vtk:
            lat.LunusWritevtk(mean_vtk_file)
        else:
            pass
        
        print "Symmetrizing lattice based on symmetry operators for Laue Class: "+self.laue_group
        lat.LunusSymlt(lunus_sym_lib[self.laue_group])
        if lat:
            lat.LunusWritelt(sym_lat_file) 
        else:
            pass
        if vtk:
            lat.LunusWritevtk(sym_vtk_file)
        else:
            pass
        
        print 'Isolating the anisotropic signal'
        lat.LunusAnisolt(self.cella, self.cellb, self.cellc, self.alpha, self.beta, self.gamma)
        if lat:
            lat.LunusWritelt(aniso_lat_file)
        else:
            pass
        if vtk:
            lat.LunusWritevtk(aniso_vtk_file)
        else:
            pass

        return

    # def average_symmetry_mates_replaced(self):

    #     ### gather files for calculations
    #     mean_vtk_file   = (work_dir+"/lattices/"+self.id+"_mean.vtk")
    #     mean_lat_file   = (work_dir+"/lattices/"+self.id+"_mean.lat")
    #     sym_lat_file    = (work_dir+"/lattices/"+self.id
    #                           +"_mean_sym.lat")
    #     sym_vtk_file    = (work_dir+"/lattices/"+self.id
    #                           +"_mean_sym.vtk")
    #     aniso_lat_file  = (work_dir+"/lattices/"+self.id+"_mean_sym_aniso.lat")
    #     aniso_vtk_file  = (work_dir+"/lattices/"+self.id+"_mean_sym_aniso.vtk")

    #     ### use Lunus methods to apply symmetry operations
    #     print "Symmetrizing lattice based on symmetry operators for Laue Class: "+self.laue_group
    #     call(['vtk2lat', mean_vtk_file, mean_lat_file])
    #     call(['symlt', mean_lat_file, sym_lat_file, str(lunus_sym_lib[self.laue_group])])
    #     print 'Isolating the anisotropic signal'
    #     call(['anisolt', sym_lat_file, aniso_lat_file, self.unit_cell])
    #     call(['lat2vtk', sym_lat_file, sym_vtk_file])
    #     call(['lat2vtk', aniso_lat_file, aniso_vtk_file])

    #     return




###--------- functional, but unnecessary --------###

def index_from_files(f_one,f_two,f_three):

    import dxtbx
    from iotbx.phil import parse
    from dxtbx.datablock import DataBlockFactory
    from dials.array_family import flex
    from dials.algorithms.indexing.indexer import indexer_base
    from dials.util.options import OptionParser
    import copy

    phil_scope_str='''
        output {{
          shoeboxes = True
        .type = bool
        .help = Save the raw pixel values inside the reflection shoeboxes.
        }}
        include scope dials.algorithms.spot_finding.factory.phil_scope
        include scope dials.algorithms.indexing.indexer.index_only_phil_scope
        include scope dials.algorithms.refinement.refiner.phil_scope
        indexing.known_symmetry.unit_cell={0}
          .type = unit_cell
        indexing.known_symmetry.space_group={1}
          .type = space_group
          '''
    phil_scope = parse(phil_scope_str.format(target_cell,target_sg), process_includes=True)
    #  from dials.util.options import OptionParser
    parser = OptionParser(phil=phil_scope)
    params, options = parser.parse_args(args=[], show_diff_phil=True)

    params.refinement.parameterisation.scan_varying = False
    params.indexing.method='real_space_grid_search'
    # params.indexing.method='fft3d'
    #  params.indexing.max_cell=800
    #  params.spotfinder.filter.min_spot_size=3
      
    filenames = [f_one,f_two,f_three]

    datablock = DataBlockFactory.from_filenames(filenames)[0]

    observed = flex.reflection_table.from_observations(datablock, params)
    observed.as_pickle("strong.pickle")
    print "Number of observed reflections:", len(observed)

    working_params = copy.deepcopy(params)
    imagesets = datablock.extract_imagesets()

    print "indexing..."
    t0 = time()
    # new dials, fix by Aaron
    idxr = indexer_base.from_parameters(observed, imagesets, params=params)
    idxr.index()
    tel = time()-t0
    print "done indexing (",tel," sec)"

    # new dials
    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments
    print experiments.crystals()[0]
    crystal_params = experiments.crystals()[0]
    with open('crystal.pkl', 'wb') as output:
        pickle.dump(crystal_params, output, pickle.HIGHEST_PROTOCOL)

    return




###----------------Direct call of script-------------------###

def main():
    exp_file, d_min = argv[1:]

    test_exp = DiffuseExperiment()
    test_exp.read_experiment_file(exp_file)

    img_set = test_exp.experiments.imagesets()
    imgs = img_set[0]
    file_list = imgs.paths()
    img_file = file_list[0]

    test_img = DiffuseImage(img_file)
    test_img.set_general_variables()
    test_img.remove_bragg_peaks(reference=True)
    test_img.scale_factor()
    test_img.crystal_geometry(test_exp.crystal)
    d_min = float(d_min)
    test_img.setup_diffuse_lattice(d_min)
    test_img.integrate_diffuse()

if __name__ == "__main__":
    main()




