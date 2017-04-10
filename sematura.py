#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

### Import required functions
import numpy as np
from subprocess import call, check_output, Popen
from time import clock, time
from glob import glob
from os.path import isdir
from sys import argv
from copy import deepcopy
import cPickle as pickle

### Import required functions specific to Phenix
from dxtbx import load
from dxtbx.format.Registry import Registry
from iotbx.phil import parse
from dxtbx.datablock import DataBlockFactory
from dials.array_family import flex
from dials.algorithms.indexing.real_space_grid_search \
     import indexer_real_space_grid_search as indexer

### Import user-defined parameters
from sematura_params import *


### Library of functions for diffuse scattering analysis
class diffuse(object):

    def __init__(self, filenumber):

        ### Define directories
        global work_dir, lunus_dir

        lunus_dir = check_output(['which', 'symlt'])
        lunus_dir = lunus_dir.replace("/c/bin/symlt\n","")
        work_dir = check_output("pwd")
        work_dir = work_dir.replace("\n","")

        ### Make subdirectories to organize output
        if (isdir(work_dir+"/proc")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/proc"])
        if (isdir(work_dir+"/radial_averages")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/radial_averages"])
        if (isdir(work_dir+"/lattices")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/lattices"])
        if (isdir(work_dir+"/arrays")) == True:
            pass
        else:
            call(['mkdir', work_dir+"/arrays"])

        ### Define filenames
        self.name   = (image_dir+"/"+image_prefix+filenumber+".img")
        self.lunus  = (work_dir+"/proc/"+lunus_image_prefix
                       +filenumber+".img")
        self.radial = (work_dir+"/radial_averages/"+lunus_image_prefix
                       +filenumber+".asc")
        self.raw    = (image_dir+"/"+image_prefix+filenumber+".cbf")
        # self.latt   = (work_dir+"/lattices/"+diffuse_lattice_prefix
        #                +"_diffuse_"+filenumber+".vtk")
        # self.counts = (work_dir+"/lattices/"+diffuse_lattice_prefix
        #                +"_counts_"+filenumber+".vtk")
        self.intensity = (work_dir+"/arrays/"+diffuse_lattice_prefix
                       +filenumber+".npz")

    ### This function converts the raw data to the img format
    def cbf2img(self, filenumber):
        
        print "Converting file #"+filenumber+" to IMG format"
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

        ### Define names of files to use for indexing
        raw_file_list = glob(image_dir+"/"+image_prefix+"*.cbf")
        indexing_data_file_one = raw_file_list[int(indexing_one)-1]
        indexing_data_file_one = indexing_data_file_one.replace(".cbf",".img")
        indexing_data_file_two = raw_file_list[int(indexing_two)-1]
        indexing_data_file_two = indexing_data_file_two.replace(".cbf",".img")
        indexing_data_file_three = raw_file_list[int(indexing_three)-1]
        indexing_data_file_three = indexing_data_file_three.replace(".cbf",".img")

        ### Input parameters for indexing
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

        ### Extract data for indexing
        datablock = DataBlockFactory.from_filenames(filenames)[0]
        observed = flex.reflection_table.from_observations(datablock, params)
        observed.as_pickle("strong.pickle")
        print "Number of observed reflections:", len(observed)
        working_params = deepcopy(params)
        imagesets = datablock.extract_imagesets()

        ### Report if beam parameters differ between images
        print "Beam parameters from image:"
        print imagesets[0].get_beam()
        print "Files used for indexing have same beam parameters:"
        print imagesets[0].get_beam() == imagesets[0].get_beam()
        print imagesets[1].get_beam() == imagesets[0].get_beam()
        print imagesets[2].get_beam() == imagesets[0].get_beam()

        ### Index the data to get the crystal matrices
        print "Indexing the data using a target space unit cell of:\n"+target_cell+"\nand a target space group number of:\n"+target_sg
        t0 = clock()
        idxr = indexer(observed, imagesets, params=working_params)
        tel = clock()-t0
        print "Done indexing (",tel," sec)"

        ### Report and store crystal matrices for mapping data to diffuse lattice
        indexed = idxr.refined_reflections
        experiments = idxr.refined_experiments
        print experiments.crystals()[0]
        crystal_params = experiments.crystals()[0]
        with open('crystal.pkl', 'wb') as output:
            pickle.dump(crystal_params, output, pickle.HIGHEST_PROTOCOL)

        return

    ### This function uses Lunus methods to remove the Bragg peaks from
    ### the raw data
    def debragg(self, filenumber):
        
        print "Removing beamstop shadow from image"
        call(['punchim', self.name, punchim_xmin, punchim_xmax,
                          punchim_ymin, punchim_ymax, filenumber+"_tmp1.img"])
        print "Cleaning up edge of image"
        call(['windim', filenumber+'_tmp1.img', windim_xmin,
                          windim_xmax, windim_ymin, windim_ymax,
                          filenumber+'_tmp2.img'])
        print "Setting detection threshold"
        call(['thrshim', filenumber+'_tmp2.img', thrshim_min,
                          thrshim_max, filenumber+'_tmp3.img'])
        print "Correcting for beam polarization"
        call(['polarim', filenumber+'_tmp3.img',
                         filenumber+'_tmp4.img', polarim_dist,
                         polarim_polarization, polarim_offset])
        print "Solid-angle normalization of image"
        call(['normim', filenumber+'_tmp4.img',
                         filenumber+'_tmp5.img', normim_tilt_x, normim_tilt_y])
        print "Removing Bragg peaks from image"
        call(['modeim', filenumber+'_tmp5.img',
                         filenumber+'_tmp6.img', modeim_kernel_width,
                         modeim_bin_size])
        print "Cleaning up directory"
        call(['cp', filenumber+'_tmp6.img', self.lunus])
        call(['rm', filenumber+"_tmp1.img", filenumber+'_tmp2.img',
                    filenumber+'_tmp3.img', filenumber+'_tmp4.img',
                    filenumber+'_tmp5.img', filenumber+'_tmp6.img'])

        return

    ### This function uses Lunus methods to create a radially averaged copy
    ### of the image for scaling purposes
    def radial_avg(self, filenumber):

        print 'Creating radially averaged copy of image #'+filenumber
        call(['avgrim', self.lunus, filenumber+'_tmp.rf'])
        infile,outfile = filenumber+'_tmp.rf',self.radial
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = Popen(
                    ['binasc', '2'],stdout=ouf,stdin=inf)
                proc.wait()

        return

    ### This function uses Lunus methods to calculate a reference statistic
    ### for calculating the scaling factor for all images
    def make_radial_ref(self):

        ### Uses the radially averaged copy of the user chosen reference
        print 'Making reference statistic for scaling'
        ref_one = open(self.radial, 'r')
        ref_one_lines = ref_one.readlines()
        ref_one.close()
        
        ### Only use ring around image defined by user
        ref_two = ref_one_lines[scale_inner_radius:scale_outer_radius]
        out_file = open('ref.asc', 'w')
        for line in ref_two:
            out_file.write(line)
        out_file.close()

        ### Create files used by scale_image calculations
        infile,outfile = 'ref.asc', 'reference.rf'
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = Popen(
                    ['binasc', '3'],stdout=ouf,stdin=inf)
                proc.wait()
        call(['mulrf', 'reference.rf', 'reference.rf', 'xx.rf'])

        return

    ### This function uses Lunus methods to calculate a scaling factor for
    ### each diffraction image which will be used during integration onto
    ### the diffuse lattice
    def scale_image(self, filenumber):

        global scale_factor, scale_factor_error

        ### Uses the radially averaged copy of the image
        print 'Calculating scaling factor for image #'+filenumber
        ref_one = open(self.radial, 'r')
        ref_one_lines = ref_one.readlines()
        ref_one.close()

        ### Only use ring around image defined by user
        ref_two = ref_one_lines[scale_inner_radius:scale_outer_radius]
        out_file = open(filenumber+'.asc', 'w')
        for line in ref_two:
            out_file.write(line)
        out_file.close()

        infile,outfile = filenumber+'.asc', filenumber+'.rf'
        with open(outfile,'w') as ouf:
            with open(infile,'r') as inf:
                proc = Popen(
                    ['binasc', '3'],stdout=ouf,stdin=inf)
                proc.wait()

        ### Use Lunus methods to find average value within ring
        call(['mulrf', 'reference.rf', filenumber+'.rf',
                         filenumber+'_xy.rf'])
        call(['mulrf', filenumber+'.rf', filenumber+'.rf',
                         filenumber+'_yy.rf'])
        xx = float(check_output(['avgrf', 'xx.rf']))
        xy = float(check_output(['avgrf', filenumber+'_xy.rf']))
        yy = float(check_output(['avgrf', filenumber+'_yy.rf']))

        ### Calculate scale factor for use during integration
        global scale_factor_error, scale_factor_error
        scale_factor = xx / xy
        scale_factor_error = np.sqrt(xx+yy*scale_factor*scale_factor
                                     -2.*scale_factor*xy)/np.sqrt(xx)
        print "This image has a scale factor of "+str(scale_factor)
        print "This image has a scale factor error of "+str(scale_factor_error)

        return

    ### This function maps diffuse data from an image onto a 3D lattice
    def procimg(self, Isize1,Isize2,scale_factor,mask_tag,A_matrix,rvec,
                DATA,latxdim,latydim,latzdim):
        
        global lat
        tmid = clock()

        ### define the lattice indices at which h,k,l = 0,0,0
        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1

        ### total number of voxels in the lattice
        latsize = latxdim*latydim*latzdim

        ### generate lattice to store data & counts for averaging
        lat = np.zeros(latsize*2, dtype=np.float32).reshape((2,latzdim,latydim,latxdim))

        ### fetch A_matrix from Dials output after indexing
        a_mat = np.array(A_matrix)
        a_mat.shape = (3,3)

        ### calculate h,k,l for every data point
        print "Calculating h,k,l for each data point from image"
        H = np.tensordot(a_mat,rvec,(1,1))
        H = np.transpose(H)
        tel = str(clock()-tmid)
        print "Calculated h,k,l for ",len(H[:]),\
              " data points ("+tel+" seconds)"

        ### fetch data from Lunus processed image using Dials methods
        val = np.array(DATA, dtype=int)

        ### rearrange order of data points to match h,k,l matrix (H)
        val.shape = (Isize2,Isize1)
        val = np.transpose(val)
        val = val.flatten()

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
        i_int += int(i0)
        j_int += int(j0)
        k_int += int(k0)

        ### make a boolean mask to be used for data & hkls
        data_mask = np.array(val, dtype=bool)
        ### mask data that are too close to bravais lattice
        data_mask[delta_i<0.25] = False
        data_mask[delta_j<0.25] = False
        data_mask[delta_k<0.25] = False
        ### mask data outside unit cell on diffuse lattice
        data_mask[i_int<0] = False
        data_mask[i_int>=latxdim] = False
        data_mask[j_int<0] = False
        data_mask[j_int>=latydim] = False
        data_mask[k_int<0] = False
        data_mask[k_int>=latzdim] = False
        ### mask data with negative intensities
        data_mask[val<=0] = False
        ### mask data marked with mask_tag by lunus software
        data_mask[val>=mask_tag] = False
        
        ### apply mask to data and associated hkls
        val = val[data_mask]
        i_int = i_int[data_mask]
        j_int = j_int[data_mask]
        k_int = k_int[data_mask]

        ### map the data onto the diffuse lattice
        np.add.at(lat[0], [k_int,j_int,i_int], val)
        diff = lat[0].flatten()
        diff_scaled = np.multiply(diff,scale_factor)
        ### keep track of the number of data points
        ### added at each lattice point (for averaging)
        val[val!=0]=1
        np.add.at(lat[1], [k_int,j_int,i_int], val)
        coun = lat[1].flatten()

        # diff_scaled[coun==0] = 32767

        lat = [diff_scaled,coun]

        return lat

    ### This function integrates the diffuse scattering using the procimg fxn
    def integrator(self):

        global lt, ct, latxdim, latydim, latzdim, i0, j0, k0, latsize

        ### set resolution of diffuse lattice
        res = float(resolution)
        latxdim = (int(cella/res)+1)*2
        latydim = (int(cellb/res)+1)*2
        latzdim = (int(cellc/res)+1)*2
        latsize = latxdim*latydim*latzdim
        print "Lattice size = ",latsize
        lt = np.zeros(latsize, dtype=np.float32)
        ct = np.zeros(latsize, dtype=np.float32)
        ### set origin and tag value for null data
        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1
        mask_tag = 32767

        ### use Dials to get data from processed image
        imgname = self.lunus
        print "Integrating file %s with scale factor %f"%(imgname,scale_factor)

        img = load(imgname)
        detector = img.get_detector()
        print detector
        beam = img.get_beam()
        scan = img.get_scan()
        gonio = img.get_goniometer()

        ### map data from pixels to mm positions
        print "mapping data from pixels to mm positions"
        t0 = clock()
        lab_coordinates = flex.vec3_double()
        for panel in detector: 
            pixels = flex.vec2_double(panel.get_image_size())
            mms = panel.pixel_to_millimeter(pixels)
            lab_coordinates.extend(panel.get_lab_coord(mms))

        s1_vectors = lab_coordinates.each_normalize() * (1/beam.get_wavelength())
        x_vectors = s1_vectors - beam.get_s0()
        print "there are ",x_vectors.size()," elements in x_vectors"

        ### get A_matrix from indexing fxn
        with open('crystal.pkl', 'rb') as input:
            crystal_params = pickle.load(input)
        crystal = deepcopy(crystal_params)
        axis = gonio.get_rotation_axis()
        start_angle, delta_angle = scan.get_oscillation()
        crystal.rotate_around_origin(axis, start_angle + (delta_angle/2), deg=True)
        A_matrix = crystal.get_A().inverse()

        ### integrate diffuse scattering from single or multi panel detectors
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
        print "done integrating diffuse scattering (",tel," sec)"

        ### divide intensities from counts
        t0 = clock()
        lt = np.add(lt,latit[0])
        ct = np.add(ct,latit[1])
        tel = clock()-t0
        print "Took ",tel," secs to update the lattice"

        np.savez(self.intensity, lt=lt, ct=ct)

        return

    ### This function writes a vtk file from a numpy array
    def array2vtk(self, in_file, key):

        array = np.load(in_file)
        array_values = array[key]
        out_file = in_file.replace("/arrays/","/lattices/")
        out_file = out_file.replace(".npz",".vtk")
        vtkfile = open(out_file,"w")

        array_values = array_values.flatten()

        res = float(resolution)
        latxdim = (int(cella/res)+1)*2
        latydim = (int(cellb/res)+1)*2
        latzdim = (int(cellc/res)+1)*2
        latsize = latxdim*latydim*latzdim
        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1
        a_recip = 1./cella
        b_recip = 1./cellb
        c_recip = 1./cellc

        print >>vtkfile,"# vtk DataFile Version 2.0"
        print >>vtkfile,"lattice_type=PR;unit_cell={0},{1},{2},{3},{4},{5};space_group={6};".format(cella,cellb,cellc,alpha,beta,gamma, target_sg)
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
                    print >>vtkfile,array_values[index],
                    index += 1
                print >>vtkfile,""

        vtkfile.close()
        return

    ### calculate mean lattice using stored numpy arrays
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
        mean_lt = np.divide(sum_int_masked,sum_ct_masked)
        mean_lt.set_fill_value(-32768)
        mean_lt_clean = mean_lt.filled()
        out = (work_dir+"/arrays/"+diffuse_lattice_prefix+"_mean.npz")
        np.savez(out, mean_lt=mean_lt_clean)

        return

    ### This function uses symmetry operations within Lunus to expand
    ### the mean lattice intensities to P1 point-group symmetry.
    def symmetrize(self):

        ### gather files for calculations
        mean_vtk_file   = (work_dir+"/lattices/"+diffuse_lattice_prefix+"_mean.vtk")
        mean_lat_file   = (work_dir+"/lattices/"+diffuse_lattice_prefix+"_mean.lat")
        sym_lat_file    = (work_dir+"/lattices/"+diffuse_lattice_prefix
                              +"_mean_sym.lat")
        sym_vtk_file    = (work_dir+"/lattices/"+diffuse_lattice_prefix
                              +"_mean_sym.vtk")
        aniso_lat_file  = (work_dir+"/lattices/"+diffuse_lattice_prefix+"_mean_sym_aniso.lat")
        aniso_vtk_file  = (work_dir+"/lattices/"+diffuse_lattice_prefix+"_mean_sym_aniso.vtk")
        ### file mapping space group to lunus input
        sg_conv = open(lunus_dir+"/analysis/sg_pg_lunus.csv","r")
        sg_conv_lines = sg_conv.readlines()
        ### parse file without pandas
        lines = None
        lunus_key = ""
        laue_class = ""
        for line in sg_conv_lines:
            line = line.split('\r')
            lines = line
        ### get laue class and Lunus input from file
        for item in lines:
            item = item.split(',')
            if item[0] == target_sg:
                lunus_key += item[4]
                laue_class += item[3]


        ### use Lunus methods to apply symmetry operations
        print "Symmetrizing lattice based on symmetry operators for Laue Class: "+laue_class
        call(['vtk2lat', mean_vtk_file, mean_lat_file])
        # call(['xflt', aniso_lat_file, 'tmp2_xf.lat', '1'])
        call(['symlt', mean_lat_file, sym_lat_file, lunus_key])
        call(['lat2vtk', sym_lat_file, sym_vtk_file])


        ### Isolate the anisotropic signal from mean lattice
        print 'Isolating the anisotropic signal'
        ### the step below decreased correlation coefficient by less than 0.001 but don't need it
        # call(['xflt', mean_lat_file, 'tmp_xf.lat', '1'])
        call(['avgrlt', sym_lat_file, 'tmp_xf.rf'])
        call(['subrflt', 'tmp_xf.rf', sym_lat_file, aniso_lat_file])
        call(['lat2vtk', aniso_lat_file, aniso_vtk_file])


        return

        



###----------------Direct call of script-------------------###

def main():
    script, filenum = argv
    filenum = float(filenum)
    filenum = '{:05.0f}'.format(filenum)

    data = diffuse(filenum)
    data.cbf2img(filenum)
    data.debragg(filenum)
    data.radial_avg(filenum)
    data.scale_image(filenum)
    # data.indexing()
    data.integrator()
    # data.latout()
    # data.ctout()

if __name__ == "__main__":
    main()




