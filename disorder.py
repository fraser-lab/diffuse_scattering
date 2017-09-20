"""
The purpose of this library is to develop automated methods
for generating models of disorder that give rise to diffuse
X-ray scattering.

The aim is to take PDB files as input and export hkl files
of electron density attributed to diffuse scattering.

Author(s):
Alexander M. Wolff      (primary author)
Michael E. Wall         (primary author)
Andrew Van Benschoten   (contributed initial methods)

"""


from __future__ import division

lunus_sym_lib = {'-1': 0, '2/m': -1, 'mmm': -2, '4/m': -3, '4/mmm': -4, \
'-3': -5, '-3m': -6, '6/m': -7, '6/mmm': -8, 'm-3': -9, 'm-3m': -10}

class DisorderModel(object):

    def __init__(self, filepath, dmin):
        
        self.filepath = filepath
        self.filename = self.filepath.split("/")[-1]
        self.name = self.filename.rstrip(".pdb")
        self.structure_factors = None
        self.hkl = None
        self.vtk = None
        self.lat = None
        self.dmin = dmin
        self.unit_cell = None
        self.space_group = None
        self.laue_group = None

        print("building a disorder model from {}".format(self.filename))


    def calculate_structure_factors(self):

        print("calculating structure factors for {}".format(self.name))
        from iotbx.pdb import hierarchy
        pdb_in = hierarchy.input(file_name=self.filepath)
        xrs = pdb_in.input.xray_structure_simple()
        fcalc = xrs.structure_factors(d_min=self.dmin).f_calc()
        fc_square = fcalc.as_intensity_array()
        fc_square_p1 = fc_square.expand_to_p1()
        self.structure_factors = fc_square_p1.generate_bijvoet_mates()
        uc = self.structure_factors.unit_cell()
        uc_nospace = str(uc).replace(" ", "")
        uc_nospace_noparen = uc_nospace[1:-1]
        self.unit_cell = uc_nospace_noparen
        self.space_group = fc_square.space_group()
        self.laue_group = self.space_group.laue_group_type()

        print("unit cell: {}".format(self.unit_cell))

        return

    def write_hkl(self):
        ### best to use internal CCTBX filewriter, but:

        ### this works but isn't the HKL we want
        #f = open("Icalc.hkl","w")
        #fc_square_p1.export_as_cns_hkl(file_name="Icalc.hkl")
        #f.close()

        ### This works, but requires scaling first, or else the file is very messy & arbitrarily scaled
        ### figure out scaling!
        #f = open("data_notcns.hkl", "w")
        #fc_square_p1_friedel.export_as_shelx_hklf(file_object=f,normalise_if_format_overflow=True)
        #f.close()

        ### So...hacky workaround:
        self.hkl = (self.name+"_Icalc.hkl")
        print("writing full set of calculated structure factors to {}".format(self.hkl))
        f = open(self.hkl,'w')
        for hkl,intensity in self.structure_factors:
            print >>f, "%4d %4d %4d    %10.2f" %(hkl+tuple((intensity,)))
        f.close()
        return


    def write_vtk(self):

        self.vtk = (self.name+"_Icalc.vtk")
        vtkfile = open(self.vtk,"w")
        ### crystal parameters
        cella, cellb, cellc, alpha, beta, gamma= self.unit_cell.parameters()
        a_recip, b_recip, c_recip = self.unit_cell.reciprocal_parameters()[:3]
        sgi = self.space_group.info()
        sg = sgi.symbol_and_number()

        ### grid dimensions for VTK file
        latxdim = (int(cella/self.dmin)+1)*2
        latydim = (int(cellb/self.dmin)+1)*2
        latzdim = (int(cellc/self.dmin)+1)*2
        latsize = latxdim*latydim*latzdim

        ### VTK cell origin
        i0=latxdim/2-1
        j0=latydim/2-1
        k0=latzdim/2-1

        ### write header for VTK file
        print >>vtkfile,"# vtk DataFile Version 2.0"
        print >>vtkfile,"lattice_type=PR;unit_cell={0},{1},{2},{3},{4},{5};space_group={6};".format(cella,cellb,cellc,alpha,beta,gamma,sg)
        print >>vtkfile,"ASCII"
        print >>vtkfile,"DATASET STRUCTURED_POINTS"
        print >>vtkfile,"DIMENSIONS %d %d %d"%(latxdim,latydim,latzdim)
        print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
        print >>vtkfile,"ORIGIN %.8f %.8f %.8f" %(-i0*a_recip,-j0*b_recip,-k0*c_recip)
        print >>vtkfile,"POINT_DATA %d"%(latsize)
        print >>vtkfile,"SCALARS volume_scalars float 1"
        print >>vtkfile,"LOOKUP_TABLE default\n"

        ### write the data to VTK file
        index = 0
        for k in range(0,latzdim):
            for j in range(0,latydim):
                for i in range(0,latxdim):
                    print >>vtkfile,self.structure_factors.data()[index],
                    index += 1
                print >>vtkfile,""

        vtkfile.close()
        return


    def write_lat(self):
        import subprocess
        if self.vtk == None:
            print("Error: cannot write lat file until vtk file is written")
        else:
            self.lat = self.vtk.replace("vtk","lat")
            subprocess.call(['vtk2lat', self.vtk, self.lat])

    def hkl2lat(self):
        self.lat = self.hkl.replace("hkl","lat")
        import subprocess
        ### improper workaround...the hkl2lat needs a template, but this should be a data file
        ### or a blank lattice of defined dimensions
        template = "data.lat"
        subprocess.call(['hkl2lat', self.hkl, self.lat, template])


    def build_llm_old(self, data_in, sigma, gamma, dmin, dmax):
        print("Creating Liquid-Like Motions (LLM) model of the diffuse scattering based on {}".format(self.name))
        # call(['vtk2lat', mean_vtk_file, mean_lat_file])
        # call(['xflt', aniso_lat_file, 'tmp2_xf.lat', '1'])
        # call(['symlt', mean_lat_file, sym_lat_file, lunus_key])
        # call(['lat2vtk', sym_lat_file, sym_vtk_file])
        # dmin = 25
        # dmax = 0
        fft_direction = {'forward': 1, 'reverse': -1}
        import subprocess
        ### make lattice of phase factors
        subprocess.call(['constlt', data_in, zero.lat, 0])
        ### calculate patterson of data
        subprocess.call(['fftlt', data_in, zero.lat, real.lat, imag.lat, fft_direction['forward']])
        ### calculate smearing lattice based on correlation length of gamma
        subprocess.call(['liquidcorrlt', data_in, corr.lat, gamma])
        ### calculate patterson of smearing lattice
        subprocess.call(['fftlt', corr.lat, zero.lat, corr_real.lat, corr_imag.lat, fft_direction['forward']])
        ### complex multiplication of patterson lattices
        subprocess.call(['mullt', real.lat, corr_real.lat, xx.lat])
        subprocess.call(['mullt', imag.lat, corr_imag.lat, yy.lat])
        subprocess.call(['sublt', xx.lat, yy.lat, mul_real.lat])
        subprocess.call(['mullt', real.lat, corr_imag.lat, xy.lat])
        subprocess.call(['mullt', imag.lat, corr_real.lat, yx.lat])
        subprocess.call(['sumlt', xy.lat, yx.lat, mul_imag.lat])
        subprocess.call(['mullt', imag.lat, corr_real.lat, yx.lat])
        ### back to reciprocal space
        subprocess.call(['fftlt', mul_real.lat, mul_imag.lat, llm.lat, tmp.lat, fft_direction['reverse']])
        ### make prefactor lattice to apply sigma, amplitude of atomic motion
        subprocess.call(['liquidfaclt', llm.lat, liquidfac.lat, sigma])
        ### apply sigma lattice to convolution product
        subprocess.call(['mullt', llm.lat, liquidfac.lat, tmp.lat])
        ### select subset of lattice based on radius of hkl
        subprocess.call(['culllt', tmp.lat, final_llm.lat, dmax, dmin])
        # ### isolate anisotropic portion of calculated lattice
        subprocess.call(['avgrlt', final_llm.lat, final_llm.rf])
        subprocess.call(['subrflt', final_llm.rf, final_llm.lat, final_llm_aniso.lat])
        # ### clean directory
        subprocess.call(['rm', tmp.lat, real.lat, imag.lat, mul_real.lat, mul_imag.lat, corr.lat,\
         corr_real.lat, corr_imag.lat, xx.lat, xy.lat, yx.lat, yy.lat, llm.lat])

        return


    def build_llm(self, x):
        import numpy as np
        print("Creating Liquid-Like Motions (LLM) model of the diffuse scattering based on {}".format(self.name))
        gamma = x[0]
        sigma = x[1]

        import subprocess

        ### make llm model
        subprocess.call(['llmlt', self.lat, self.name+"_llm.lat", self.unit_cell, str(gamma), str(sigma)])
        subprocess.call(['symlt', self.name+"_llm.lat", self.name+"_llm_sym.lat", str(lunus_sym_lib[self.laue_group])])
        subprocess.call(['anisolt', self.name+"_llm_sym.lat", self.name+"_llm_sym_aniso.lat", self.unit_cell])
        subprocess.call(['cullreslt', self.name+"_llm_sym_aniso.lat", self.name+"_llm_sym_aniso_culled.lat", "31.2", "1.45", self.unit_cell])
        R = subprocess.check_output(['rfaclt', "data_aniso_culled.lat", self.name+"_llm_sym_aniso_culled.lat"])
        R = R[3:10]
        print("R factor of current LLM model = {}".format(R))
        R = float(R)

        CC = subprocess.check_output(['corrlt', "data_aniso_culled.lat", self.name+"_llm_sym_aniso_culled.lat"])
        # R = R[3:10]
        print("CC factor of current LLM model = {}".format(CC))
        CC = -1*float(CC)

        return R


    def build_rigid_body_translation(self, x):
        print("Creating Rigid Body Translation (RBT) model of the diffuse scattering based on {}".format(self.name))

        sigma = x

        import subprocess
        ### make rbt model
        subprocess.call(['rbtlt', self.lat, self.name+"_rbt.lat", self.unit_cell, str(sigma)])
        subprocess.call(['symlt', self.name+"_rbt.lat", self.name+"_rbt_sym.lat", str(lunus_sym_lib[self.laue_group])])
        subprocess.call(['anisolt', self.name+"_rbt_sym.lat", self.name+"_rbt_sym_aniso.lat", self.unit_cell])
        subprocess.call(['cullreslt', self.name+"_rbt_sym_aniso.lat", self.name+"_rbt_sym_aniso_culled.lat", "31.2", "1.45", self.unit_cell])

        R = subprocess.check_output(['rfaclt', "data_aniso_culled.lat", self.name+"_rbt_sym_aniso_culled.lat"])
        R = R[3:10]
        print("R factor of current RBT model = {}".format(R))
        R = float(R)

        CC = subprocess.check_output(['corrlt', "data_aniso_culled.lat", self.name+"_rbt_sym_aniso_culled.lat"])
        # R = R[3:10]
        print("CC factor of current RBT model = {}".format(CC))
        CC = -1*float(CC)

        return R


    def build_nma(self):
        print("NMA model invoked")
        from prody import *
        from pylab import *
        cypa = parsePDB('5f66')
        calphas = cypa.select('protein and name CA')
        anm = ANM('cypa ANM analysis')
        anm.buildHessian(calphas, cutoff=25, gamma=10.5)
        print("cutoff:{}".format(anm.getCutoff()))
        print("gamma:{}".format(anm.getGamma()))
        anm.calcModes(10)
        cov = anm.getCovariance().round(2)
        import matplotlib.pyplot as plt

        # print(cov)
        plt.imshow(cov)
        plt.show()

        return

    def refined_llm(self):

        import scipy.optimize as optimize
        rranges = ((1, 20), (0.1, 2.0))
        # rranges = (7.5,0.5)

        # best_params = optimize.differential_evolution(self.build_llm, rranges)
        output = optimize.brute(self.build_llm, rranges, Ns=20, full_output=True, finish=optimize.fmin_powell)
        # best_params = optimize.minimize(self.build_llm, rranges, method='Powell')
        # best_params = optimize.basinhopping(self.build_llm, rranges)
        # gamma = best_params[0][0]
        # sigma = best_params[0][1]
        # R_factor = best_params[1]
        # print("\nRefined gamma = {}".format(gamma))
        # print("\nRefined sigma = {}".format(sigma))
        # print("\nR factor for refined model = {}".format(R_factor))
        best_params = output[0]
        best_score  = output[1]
        param_space = output[2]
        param_space_scores = output[3]

        print("best parameters:\n{}".format(best_params))
        print("minimized R-factor:\n{}".format(best_score))
        import matplotlib.pyplot as plt
        fig = plt.figure(num=1,figsize=(15,15))
        CS = plt.contour(param_space[0],param_space[1],param_space_scores, cmap='plasma')
        plt.clabel(CS, inline=1, fontsize=20)
        plt.xlabel("$\gamma$ (ang)", fontsize=20)
        plt.ylabel("$\sigma$ (ang)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.title("Finding the Optimal LLM Model for 5f66: R Factor Distribution", fontsize=30)
        plt.tight_layout()
        plt.savefig("optimized_llm_5f66.png")


        return


    def refined_rbt(self):

        import scipy.optimize as optimize
        rranges = [(0.0,5.0)]

        output = optimize.brute(self.build_rigid_body_translation, rranges, full_output=True, finish=optimize.fmin)
        # best_params = optimize.differential_evolution(self.build_rigid_body_translation, rranges)
        # best_params = optimize.minimize(self.build_rigid_body_translation, rranges, method='Powell')
        # sigma = best_params[0]
        # R_factor = best_params[1]
        # print("\nRefined sigma = {}".format(sigma))
        # print("\nR factor for refined model = {}".format(R_factor))
        best_params = output[0]
        best_score  = output[1]
        param_space = output[2]
        param_space_scores = output[3]

        print("best parameters:\n{}".format(best_params))
        print("minimized R-factor:\n{}".format(best_score))
        # print(param_space)
        # print(param_space_scores)
        import matplotlib.pyplot as plt
        fig = plt.figure(num=1,figsize=(15,15))
        # CS = plt.contour(param_space[0],param_space[1],param_space_scores, cmap='plasma')
        # plt.clabel(CS, inline=1, fontsize=20)
        plt.plot(param_space,param_space_scores)
        plt.ylabel("R factor", fontsize=20)
        plt.xlabel("$\sigma$ (ang)", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.title("Finding the Optimal RBT Model for 5f66: R Factor Distribution", fontsize=30)
        plt.tight_layout()
        plt.savefig("optimized_rbt_5f66.png")
        print best_params

        return




### test the class

A = DisorderModel("./5f66.pdb", 1.4)


A.calculate_structure_factors()
A.write_hkl()
A.hkl2lat()
# A.write_vtk()
# A.write_lat()
# import numpy as np
# x = np.array([7.5,0.5])
# A.build_llm(x)
# A.refined_llm()
# A.build_rigid_body_translation(x)
# A.refined_rbt()
A.build_nma()
