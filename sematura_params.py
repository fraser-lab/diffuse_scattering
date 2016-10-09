
### User input for diffuse scattering analysis


###---------Parameters from Bragg data reduction---------###
alpha       =   90.0                # alpha angle of unit cell
beta        =   90.0                # beta angle of unit cell
gamma       =   90.0                # gamma angle of unit cell
cella       =   43.16               # length of unit cell along a axis
cellb       =   52.13               # length of unit cell along b axis
cellc       =   91.67               # length of unit cell along c axis
target_sg   =   str(19)             # space group from Bragg data reduction
resolution  =   str(1.4)            # maximum resolution for diffuse data integration
target_cell =   (str(cella)+" "+str(cellb)+" "+str(cellc)+" "
                +str(alpha)+" "+str(beta)+" "+str(gamma))


###---------------file and directory names---------------###
diffuse_lattice_prefix  =   'cypa'
image_prefix            =   'set_1_1_'
lunus_image_prefix      =   'set_1_1lunus_'
image_dir               =   '/netapp/home/awolff/cypa_redo/data'
phenix_dir              =   '/netapp/home/awolff/packages/phenix-dev-2499'

###-----------Choose files to use for indexing-----------###
indexing_one   = str(10)            # number of first file for indexing
indexing_two   = str(50)            # number of second file for indexing
indexing_three = str(90)            # number of third file for indexing



###----------------Lunus input parameters----------------###
### image borders
windim_xmax     =str(2363)          # right border for processed image (pixels)
windim_xmin     =str(10)            # left border for processed image (pixels)
windim_ymax     =str(2427)          # top border for processed image (pixels)
windim_ymin     =str(100)           # bottom border for processed image (pixels)
### beamstop borders
punchim_xmax    =str(2459)          # right border of beam stop shadow (pixels)
punchim_xmin    =str(1156)          # left border of beam stop shadow (pixels)
punchim_ymax    =str(1363)          # top border of beam stop shadow (pixels)
punchim_ymin    =str(1224)          # bottom border of beam stop shadow (pixels)
### intensity thresholds
thrshim_max     =str(10000)         # maximum reliable intensity for detector (count)
thrshim_min     =str(1)             # minimum intensity for diffuse analysis (count)
### binning to remove Bragg peaks
modeim_bin_size     =str(1)         # mask bin for mode filtering (count)
modeim_kernel_width =str(20)        # mask width for mode filtering (pixels)
### beam polarization correction
polarim_dist        =str(183.02)    # distance from sample to detector (mm)
polarim_offset      =str(0.0)       # ask Mike W.
polarim_polarization=str(0.98)      # polarization fraction (0.0-1.0) 
### solid-angle image normalization
normim_tilt_x       =str(0.0)       # ask Mike W.
normim_tilt_y       =str(0.0)       # ask Mike W.
### image scaling
reference_image_number =str(1)     # reference for scaling
scale_inner_radius     =100        # ask Mike W.
scale_outer_radius     =800        # ask Mike W.


