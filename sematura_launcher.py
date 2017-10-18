#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

### Import required functions
from subprocess import check_output, call
from os import getcwd
from os.path import dirname
from glob import glob
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sys import exit, argv
from sematura import DiffuseExperiment, DiffuseImage

# ### Import user-defined parameters & log them
# from sematura_params import *

# def logger():

#     import sematura_params
    
#     paramlog = open('param.log', 'w')

#     for item in dir(sematura_params):
#         if item[0] != '_':
#             line = str(item) + ' = ' + str(getattr(sematura_params, item)) + '\n'
#             paramlog.write(line)

#     return



### Setup command-line input flags and help messages
epi ='''
Dependencies:
    User must have DIALS software package installed & sourced.

Control flow:
    Use "i", "p", and "a" flags. All three can be launched in a single job.
    eg: libtbx.python sematura_launcher.py -i -p -a
    '''

parser = ArgumentParser(usage='libtbx.python sematura_launcher.py [options]', description='Analyze diffuse scattering in diffraction images.', epilog=epi, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i', '--init', help='calculate reference statistic & prepare indexing files (must perform before or with processing step)', action='store_true')
parser.add_argument('-p', '--process', help='processes all images and maps diffuse scattering to reciprocal lattice', action='store_true')
parser.add_argument('-a', '--analysis', help='calculates mean & symmetrized mean lattices', action='store_true')
parser.add_argument('-n', '--nocluster', help='allows user to run on non SGE environment (default submits jobs to SGE)', action='store_true')
parser.add_argument('-exp', '--experiment_file', help='input experiment file from DIALS or XDS')
parser.add_argument('-d', '--d_min', help='input highest resolution bin for data processing')
parser.add_argument('-ref', '--reflection_file', help='input reflection file from DIALS or XDS')

args = parser.parse_args()


### Variable definitions
# lunus_dir = check_output(['which', 'symlt'])
# lunus_dir = lunus_dir.replace("/c/bin/symlt\n","")
dials = check_output(["which", "dials.version"])
dials_dir = dirname(dials)
dials_setup = dials_dir.replace("bin","setpaths.sh")
work_dir = getcwd()



if not args.experiment_file:
    print "error: Must provide at experiment file from DIALS or XDS"
    parser.print_help()
    exit(1)

if not args.d_min:
    print "error: Must provide d_min"
    parser.print_help()
    exit(1)

### Remind user to set flags
if not any([item for item in [args.init, args.process, args.analysis]]):
    print "error: Must provide at least one instructional argument (-i -p -a)\n"
    parser.print_help()
    exit(1)




prep = DiffuseExperiment()
prep.read_experiment_file(args.experiment_file)
img_set = prep.experiments.imagesets()
imgs = img_set[0]
file_list = imgs.paths()
img_one =  file_list[0]
num_files = len(file_list)
tasknames = "0"
for item in file_list:
    tasknames += " " + item


### Initialize & prepare for data analysis
if args.init:
    ### Prepare reference
    initializer = DiffuseImage(img_one)
    initializer.set_general_variables()
    initializer.remove_bragg_peaks(reference=True)
    #initializer.radial_average(reference=True)


    # logger()
    # initializer = diffuse(reference_number)
    # initializer.cbf2img(reference_number)
    # initializer.debragg(reference_number)
    # initializer.radial_avg(reference_number)
    # initializer.make_radial_ref()
    # ### Prepare indexing files
    # indexer_one = diffuse(indexing_number_one)
    # indexer_one.cbf2img(indexing_number_one)
    # indexer_two = diffuse(indexing_number_two)
    # indexer_two.cbf2img(indexing_number_two)
    # indexer_three = diffuse(indexing_number_three)
    # indexer_three.cbf2img(indexing_number_three)
    # initializer.indexing()

### Process images, from raw data to diffuse lattice file for each
if args.process:
    if not args.nocluster:

        ### Make an SGE submission script
        f = open("run_sematura.sh", "w")
        f.write("""
#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l mem_free=1G
#$ -l h_rt=2:00:00
#$ -t 1-%d
hostname
date
tasks=(%s)
input="${tasks[$SGE_TASK_ID]}"

cd %s

. %s

libtbx.python %s/scripts/sematura.py %s $input %s

        """ % (work_dir, work_dir, num_files, tasknames, work_dir, dials_setup, work_dir, args.experiment_file, args.d_min))
        f.close()

        ### submit jobs to SGE cluster
        call(['qsub', '-sync', 'y', 'run_sematura.sh'])

### Alternative processing for use without cluster environment
if args.process & args.nocluster:
    for item in file_list:
        call(['libtbx.python', work_dir+'/scripts/sematura.py', args.experiment_file, item, args.d_min])

### Analyze diffuse lattices and expand to P1 symmetry
if args.analysis:
    end = DiffuseExperiment()
    end.read_experiment_file(args.experiment_file)
    analyzer = DiffuseImage(img_one)
    analyzer.set_general_variables()
    analyzer.crystal_geometry(end.crystal)
    analyzer.setup_diffuse_lattice(float(args.d_min))
    analyzer.mean_lattice()
    analyzer.array_to_vtk()
    analyzer.average_symmetry_mates(vtk=True, lat=True)


