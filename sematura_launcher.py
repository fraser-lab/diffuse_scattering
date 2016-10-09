#---------------------------------------------------------------
###  =  Notes
#    =  code that has been commented out
#---------------------------------------------------------------

### Import required modules
from subprocess import check_output, call
from os import getcwd
from glob import glob
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sys import exit, argv
from sematura import diffuse

### Import user-defined parameters
from sematura_params import *


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
# print(args.accumulate(args.integers))



###----------------Begin Variable Definitions-------------------###

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

###-----------------End Variable Definitions--------------------###


###---------Runs fxns within Diffuse class...all work------------###
if len(argv)==1:
    print "Error: Must provide at least one argument listed below\n"
    parser.print_help()
    exit(1)

if args.init:
    ### Prepare reference
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
    # reference.indexing()


if args.process:
    # Make an SGE submission file for all of them
      
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


    call(['qsub', '-sync', 'y', 'run_sematura.sh'])

if args.process & args.nocluster:
    for item in files:
        call(['libtbx.python', lunus_dir+'/scripts/sematura.py', item])

if args.analysis:
    analyzer = diffuse(reference_image_number)
    analyzer.mean_lattice()
    analyzer.symmetrize()


