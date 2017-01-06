#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

### Import required functions
from subprocess import check_output, call
from os import getcwd
from glob import glob
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sys import exit, argv
from sematura import diffuse

### Import user-defined parameters & log them
from sematura_params import *

def logger():

    import sematura_params
    
    paramlog = open('param.log', 'w')

    for item in dir(sematura_params):
        if item[0] != '_':
            line = str(item) + ' = ' + str(getattr(sematura_params, item)) + '\n'
            paramlog.write(line)

    return



### Setup command-line input flags and help messages
epi ='''
Dependencies:
    User must have LUNUS software package installed and compiled
    User must have PHENIX software package installed & sourced.

Control flow:
    Use "i", "p", and "a" flags. All three can be launched in a single job.
    eg: libtbx.python sematura_launcher.py -i -p -a
    '''

parser = ArgumentParser(usage='libtbx.python sematura_launcher.py [options]', description='Analyze diffuse scattering in diffraction images.', epilog=epi, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i', '--init', help='calculate reference statistic & prepare indexing files (must perform before or with processing step)', action='store_true')
parser.add_argument('-p', '--process', help='processes all images and maps diffuse scattering to reciprocal lattice', action='store_true')
parser.add_argument('-a', '--analysis', help='calculates mean & symmetrized mean lattices', action='store_true')
parser.add_argument('-n', '--nocluster', help='allows user to run on non SGE environment (default submits jobs to SGE)', action='store_true')
args = parser.parse_args()


### Variable definitions
lunus_dir = check_output(['which', 'symlt'])
lunus_dir = lunus_dir.replace("/c/bin/symlt\n","")
work_dir = check_output("pwd")
work_dir = work_dir.replace("\n","")
raw_file_list = glob(image_dir+"/"+image_prefix+"*.cbf")
num_files = len(raw_file_list)
tasknames = "0"
for item in raw_file_list:
    f, t = item.split(".")
    tasknames += " " + f[-5:]

files = tasknames.split(' ')
reference_number      = files[int(reference_image_number)]
indexing_number_one   = files[int(indexing_one)]
indexing_number_two   = files[int(indexing_two)]
indexing_number_three = files[int(indexing_three)]


### Remind user to set flags
if len(argv)==1:
    print "error: Must provide at least one argument listed below\n"
    parser.print_help()
    exit(1)

### Initialize & prepare for data analysis
if args.init:
    ### Prepare reference
    logger()
    initializer = diffuse(reference_number)
    initializer.cbf2img(reference_number)
    initializer.debragg(reference_number)
    initializer.radial_avg(reference_number)
    initializer.make_radial_ref()
    ### Prepare indexing files
    indexer_one = diffuse(indexing_number_one)
    indexer_one.cbf2img(indexing_number_one)
    indexer_two = diffuse(indexing_number_two)
    indexer_two.cbf2img(indexing_number_two)
    indexer_three = diffuse(indexing_number_three)
    indexer_three.cbf2img(indexing_number_three)
    initializer.indexing()

### Process images, from raw data to diffuse lattice file for each
if args.process:

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

. %s/phenix_env.sh

libtbx.python %s/scripts/sematura.py $input

    """ % (getcwd(), getcwd(), num_files, tasknames, work_dir, phenix_dir, lunus_dir))
    f.close()

    ### submit jobs to SGE cluster
    call(['qsub', '-sync', 'y', 'run_sematura.sh'])

### Alternative processing for use without cluster environment
if args.process & args.nocluster:
    for item in files:
        call(['libtbx.python', lunus_dir+'/scripts/sematura.py', item])

### Analyze diffuse lattices and expand to P1 symmetry
if args.analysis:
    analyzer = diffuse(reference_image_number)
    analyzer.mean_lattice()
    analyzer.array2vtk(work_dir+"/arrays/"+diffuse_lattice_prefix+"_mean.npz", "mean_lt")
    analyzer.symmetrize()


