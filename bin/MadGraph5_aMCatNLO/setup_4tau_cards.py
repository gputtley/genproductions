# Phenomonoly closely follows https://arxiv.org/pdf/2104.10175.pdf

import math as m
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--generate', help= 'Generate gridpacks',  action='store_true')
parser.add_argument('--submit', help= 'Submit gridpack generation jobs',  action='store_true')
args = parser.parse_args()


#mA_list = [100,300,500]
#mphi_list = [100,300,500]
mA_list = [500]
mphi_list = [500]


### Input params needed for type X 2HDM ###
tanb   = 1.00
vev    = 246.0
yme    = 0.00
ymm    = 0.00
ymtau  = 1.777
ymup   = 0.00
ymdo   = 0.00
ymc    = 0.00
yms    = 0.00
ymt    = 173.00
ymb    = 4.70

### Calculated values ###

GLR1x1 = (yme/vev) * (2**0.5) * tanb # Type I use -/, Type II use *, Type X use *
GLR2x2 = (ymm/vev) * (2**0.5) * tanb
GLR3x3 = (ymtau/vev) * (2**0.5) * tanb

GLI1x1 = 0
GLI2x2 = 0
GLI3x3 = 0

GDR1x1 = -(ymdo/vev) * (2**0.5) / tanb # Type I use -/, Type II use *, Type X use -/
GDR2x2 = -(yms/vev) * (2**0.5) / tanb
GDR3x3 = -(ymb/vev) * (2**0.5) / tanb

GDI1x1 = 0
GDI2x2 = 0
GDI3x3 = 0

GUR1x1 = -(yme/vev) * (2**0.5) / tanb # Type I use -/, Type II use -/, Type X use -/
GUR2x2 = -(ymm/vev) * (2**0.5) / tanb
GUR3x3 = -(ymtau/vev) * (2**0.5) / tanb

GUI1x1 = 0
GUI2x2 = 0
GUI3x3 = 0

l2 = (125/vev)**2 # lambda 2 SM only
l3 = (125/vev)**2 # lambda 3 SM only

### Cards for Madgraph ###

customize_card = [
"set param_card mass 25 __mh__", # mh1 - light neutral scalar Higgs mass
"set param_card mass 35 __mH__", # mh2 - heavy neutral scalar Higgs mass
"set param_card mass 36 __mA__", # mh3 - pseudoscalar Higgs mass
"set param_card mass 37 __mHc__", # mhc - charged Higgs mass
"set param_card higgs 6 0", # mixh2 - CP conserving
"set param_card higgs 7 0", # mixh3 - CP conserving
"set param_card higgs 3 0", # LR7 - CP conserving
"set param_card higgs 4 0", # LI7 - CP conserving
"set param_card higgs 5 __mixh__", # mixh - alignment scenario sin(b-a)=1 => mixh=0, cos(b-a)=1 => mixh=pi/2
"set param_card higgs 1 {}".format(l2), # lambda 2 Higgs potential, assumed SM
"set param_card higgs 2 {}".format(l3), # lambda 3 Higgs potential, assumed SM

"set param_card yukawaglr 1 1 {}".format(GLR1x1),
"set param_card yukawaglr 1 2 0",
"set param_card yukawaglr 1 3 0",
"set param_card yukawaglr 2 1 0",
"set param_card yukawaglr 2 2 {}".format(GLR2x2),
"set param_card yukawaglr 2 3 0",
"set param_card yukawaglr 3 1 0",
"set param_card yukawaglr 3 2 0",
"set param_card yukawaglr 3 3 {}".format(GLR3x3),


"set param_card yukawagli 1 1 {}".format(GLI1x1),
"set param_card yukawagli 1 2 0",
"set param_card yukawagli 1 3 0",
"set param_card yukawagli 2 1 0",
"set param_card yukawagli 2 2 {}".format(GLI2x2),
"set param_card yukawagli 2 3 0",
"set param_card yukawagli 3 1 0",
"set param_card yukawagli 3 2 0",
"set param_card yukawagli 3 3 {}".format(GLI3x3),

"set param_card yukawagdr 1 1 {}".format(GDR1x1),
"set param_card yukawagdr 1 2 0",
"set param_card yukawagdr 1 3 0",
"set param_card yukawagdr 2 1 0",
"set param_card yukawagdr 2 2 {}".format(GDR2x2),
"set param_card yukawagdr 2 3 0",
"set param_card yukawagdr 3 1 0",
"set param_card yukawagdr 3 2 0",
"set param_card yukawagdr 3 3 {}".format(GDR3x3),


"set param_card yukawagdi 1 1 {}".format(GDI1x1),
"set param_card yukawagdi 1 2 0",
"set param_card yukawagdi 1 3 0",
"set param_card yukawagdi 2 1 0",
"set param_card yukawagdi 2 2 {}".format(GDI2x2),
"set param_card yukawagdi 2 3 0",
"set param_card yukawagdi 3 1 0",
"set param_card yukawagdi 3 2 0",
"set param_card yukawagdi 3 3 {}".format(GDI3x3),

"set param_card yukawagur 1 1 {}".format(GUR1x1),
"set param_card yukawagur 1 2 0",
"set param_card yukawagur 1 3 0",
"set param_card yukawagur 2 1 0",
"set param_card yukawagur 2 2 {}".format(GUR2x2),
"set param_card yukawagur 2 3 0",
"set param_card yukawagur 3 1 0",
"set param_card yukawagur 3 2 0",
"set param_card yukawagur 3 3 {}".format(GUR3x3),


"set param_card yukawagui 1 1 {}".format(GUI1x1),
"set param_card yukawagui 1 2 0",
"set param_card yukawagui 1 3 0",
"set param_card yukawagui 2 1 0",
"set param_card yukawagui 2 2 {}".format(GUI2x2),
"set param_card yukawagui 2 3 0",
"set param_card yukawagui 3 1 0",
"set param_card yukawagui 3 2 0",
"set param_card yukawagui 3 3 {}".format(GUI3x3),
                 ]


proc_card = [
"set default_unset_couplings 99",
"set group_subprocesses Auto",
"set ignore_six_quark_processes False",
"set loop_optimized_output True",
"set loop_color_flows False",
"set gauge unitary",
"set complex_mass_scheme False",
"set max_npoint_for_channel 0",
"import model sm",
"define p = g u c d s u~ c~ d~ s~",
"define j = g u c d s u~ c~ d~ s~",
"define l+ = e+ mu+",
"define l- = e- mu-",
"define vl = ve vm vt",
"define vl~ = ve~ vm~ vt~",
"import model 2HDM_NLO",
"define p = g u c d s b u~ c~ d~ s~ b~",
"define j = p",
"generate p p > __phi__ h3 [QCD]",
"output phi__mphi__A__mA__To4Tau -nojpeg ",
            ]

run_card = [
"#***********************************************************************",
"#                        MadGraph5_aMC@NLO                             *",
"#                                                                      *",
"#                      run_card.dat aMC@NLO                            *",
"#                                                                      *",
"#  This file is used to set the parameters of the run.                 *",
"#                                                                      *",
"#  Some notation/conventions:                                          *",
"#                                                                      *",
"#   Lines starting with a hash (#) are info or comments                *",
"#                                                                      *",
"#   mind the format:   value    = variable     ! comment               *",
"#                                                                      *",
"#   Some of the values of variables can be list. These can either be   *",
"#   comma or space separated.                                          *",
"#                                                                      *",
"#   To display additional parameter, you can use the command:          *",
"#      update to_full                                                  *",
"#***********************************************************************",
"#",
"#*******************                                                 ",
"# Running parameters",
"#*******************                                                 ",
"#",
"#***********************************************************************",
"# Tag name for the run (one word)                                      *",
"#***********************************************************************",
"  tag_1     = run_tag ! name of the run",
"#***********************************************************************",
"# Number of LHE events (and their normalization) and the required      *",
"# (relative) accuracy on the Xsec.                                     *",
"# These values are ignored for fixed order runs                        *",
"#***********************************************************************",
" 10000 = nevents ! Number of unweighted events requested",
" -1.0 = req_acc ! Required accuracy (-1=auto determined from nevents)",
" -1 = nevt_job! Max number of events per job in event generation.",
"                 !  (-1= no split).",
"#***********************************************************************",
"# Normalize the weights of LHE events such that they sum or average to *",
"# the total cross section                                              *",
"#***********************************************************************",
" average = event_norm    ! valid settings: average, sum, bias",
"#***********************************************************************",
"# Number of points per itegration channel (ignored for aMC@NLO runs)   *",
"#***********************************************************************",
" 0.01   = req_acc_FO       ! Required accuracy (-1=ignored, and use the",
"                           ! number of points and iter. below)",
"# These numbers are ignored except if req_acc_FO is equal to -1",
" 5000   = npoints_FO_grid  ! number of points to setup grids",
" 4      = niters_FO_grid   ! number of iter. to setup grids",
" 10000  = npoints_FO       ! number of points to compute Xsec",
" 6      = niters_FO        ! number of iter. to compute Xsec",
"#***********************************************************************",
"# Random number seed                                                   *",
"#***********************************************************************",
" 0    = iseed       ! rnd seed (0=assigned automatically=default))",
"#***********************************************************************",
"# Collider type and energy                                             *",
"#***********************************************************************",
" 1   = lpp1    ! beam 1 type (0 = no PDF)",
" 1   = lpp2    ! beam 2 type (0 = no PDF)",
" 6500.0   = ebeam1  ! beam 1 energy in GeV",
" 6500.0   = ebeam2  ! beam 2 energy in GeV",
"#***********************************************************************",
"# PDF choice: this automatically fixes also alpha_s(MZ) and its evol.  *",
"#***********************************************************************",
" nn23nlo = pdlabel ! PDF set",
" 244600  = lhaid   ! If pdlabel=lhapdf, this is the lhapdf number. Only",
"              ! numbers for central PDF sets are allowed. Can be a list;",
"              ! PDF sets beyond the first are included via reweighting.",
"#***********************************************************************",
"# Include the NLO Monte Carlo subtr. terms for the following parton    *",
"# shower (HERWIG6 | HERWIGPP | PYTHIA6Q | PYTHIA6PT | PYTHIA8)         *",
"# WARNING: PYTHIA6PT works only for processes without FSR!!!!          *",
"#***********************************************************************",
"  HERWIG6   = parton_shower",
"  1.0 = shower_scale_factor ! multiply default shower starting",
"                                  ! scale by this factor",
"#***********************************************************************",
"# Renormalization and factorization scales                             *",
"# (Default functional form for the non-fixed scales is the sum of      *",
"# the transverse masses divided by two of all final state particles    * ",
"# and partons. This can be changed in SubProcesses/set_scales.f or via *",
"# dynamical_scale_choice option)                                       *",
"#***********************************************************************",
" False    = fixed_ren_scale  ! if .true. use fixed ren scale",
" False    = fixed_fac_scale  ! if .true. use fixed fac scale",
" 91.118   = muR_ref_fixed    ! fixed ren reference scale",
" 91.118   = muF_ref_fixed    ! fixed fact reference scale",
" -1 = dynamical_scale_choice ! Choose one (or more) of the predefined",
"           ! dynamical choices. Can be a list; scale choices beyond the",
"           ! first are included via reweighting",
" 1.0  = muR_over_ref  ! ratio of current muR over reference muR",
" 1.0  = muF_over_ref  ! ratio of current muF over reference muF",
"#*********************************************************************** ",
"# Reweight variables for scale dependence and PDF uncertainty          *",
"#***********************************************************************",
" 1.0, 2.0, 0.5 = rw_rscale ! muR factors to be included by reweighting",
" 1.0, 2.0, 0.5 = rw_fscale ! muF factors to be included by reweighting",
" True = reweight_scale ! Reweight to get scale variation using the",
"            ! rw_rscale and rw_fscale factors. Should be a list of",
"            ! booleans of equal length to dynamical_scale_choice to",
"            ! specify for which choice to include scale dependence.",
" False = reweight_PDF  ! Reweight to get PDF uncertainty. Should be a",
"            ! list booleans of equal length to lhaid to specify for",
"            !  which PDF set to include the uncertainties.",
"#***********************************************************************",
"# Store reweight information in the LHE file for off-line model-       *",
"# parameter reweighting at NLO+PS accuracy                             *",
"#***********************************************************************",
" False = store_rwgt_info ! Store info for reweighting in LHE file",
"#***********************************************************************",
"# ickkw parameter:                                                     *",
"#   0: No merging                                                      *",
"#   3: FxFx Merging - WARNING! Applies merging only at the hard-event  *",
"#      level. After showering an MLM-type merging should be applied as *",
"#      well. See http://amcatnlo.cern.ch/FxFx_merging.htm for details. *",
"#   4: UNLOPS merging (with pythia8 only). No interface from within    *",
"#      MG5_aMC available, but available in Pythia8.                    *",
"#  -1: NNLL+NLO jet-veto computation. See arxiv:1412.8408 [hep-ph].    *",
"#***********************************************************************",
" 0        = ickkw",
"#***********************************************************************",
"#",
"#***********************************************************************",
"# BW cutoff (M+/-bwcutoff*Gamma). Determines which resonances are      *",
"# written in the LHE event file                                        *",
"#***********************************************************************",
" 15.0  = bwcutoff",
"#***********************************************************************",
"# Cuts on the jets. Jet clustering is performed by FastJet.            *",
"#  - When matching to a parton shower, these generation cuts should be *",
"#    considerably softer than the analysis cuts.                       *",
"#  - More specific cuts can be specified in SubProcesses/cuts.f        *",
"#***********************************************************************",
"  1.0  = jetalgo   ! FastJet jet algorithm (1=kT, 0=C/A, -1=anti-kT)",
"  0.7  = jetradius ! The radius parameter for the jet algorithm",
" 10.0  = ptj       ! Min jet transverse momentum",
" -1.0  = etaj      ! Max jet abs(pseudo-rap) (a value .lt.0 means no cut)",
"#***********************************************************************",
"# Cuts on the charged leptons (e+, e-, mu+, mu-, tau+ and tau-)        *",
"# More specific cuts can be specified in SubProcesses/cuts.f           *",
"#***********************************************************************",
"  0.0  = ptl     ! Min lepton transverse momentum",
" -1.0  = etal    ! Max lepton abs(pseudo-rap) (a value .lt.0 means no cut)",
"  0.0  = drll    ! Min distance between opposite sign lepton pairs",
"  0.0  = drll_sf ! Min distance between opp. sign same-flavor lepton pairs",
"  0.0  = mll     ! Min inv. mass of all opposite sign lepton pairs",
" 30.0  = mll_sf  ! Min inv. mass of all opp. sign same-flavor lepton pairs",
"#***********************************************************************",
"# Photon-isolation cuts, according to hep-ph/9801442. When ptgmin=0,   *",
"# all the other parameters are ignored.                                *",
"# More specific cuts can be specified in SubProcesses/cuts.f           *",
"#***********************************************************************",
" 20.0  = ptgmin    ! Min photon transverse momentum",
" -1.0  = etagamma  ! Max photon abs(pseudo-rap)",
"  0.4  = R0gamma   ! Radius of isolation code",
"  1.0  = xn        ! n parameter of eq.(3.4) in hep-ph/9801442",
"  1.0  = epsgamma  ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442",
" True  = isoEM  ! isolate photons from EM energy (photons and leptons)",
"#***********************************************************************",
"# Cuts associated to MASSIVE particles identified by their PDG codes.  *",
"# All cuts are applied to both particles and anti-particles, so use    *",
"# POSITIVE PDG CODES only. Example of the syntax is {6 : 100} or       *",
"# {6:100, 25:200} for multiple particles                               *",
"#***********************************************************************",
"  {} = pt_min_pdg ! Min pT for a massive particle",
"  {} = pt_max_pdg ! Max pT for a massive particle",
"  {} = mxx_min_pdg ! inv. mass for any pair of (anti)particles",
"#***********************************************************************",
"# For aMCfast+APPLGRID use in PDF fitting (http://amcfast.hepforge.org)*",
"#***********************************************************************",
" 0 = iappl ! aMCfast switch (0=OFF, 1=prepare grids, 2=fill grids)",
"#***********************************************************************",
]

submit_file = [
"source /vols/grid/cms/setup.sh",
"./gridpack_generation.sh __filename__ __directory__/__filename__"
]


def AppendLinesInList(lst,search,replace):
  if type(search)==str:
    search = [search]
    replace = [replace]
  new_list = []
  for i in lst:
    outline = i[:]
    for ind,it in enumerate(search):
      outline = outline.replace(str(it),str(replace[ind]))
    new_list.append(outline)
  return new_list

def WriteListToFile(lst,output_name):
  textfile = open(output_name, "w")
  for i in lst:
    textfile.write(i + "\n")
  textfile.close()

directory = "cards/4tau"
if not os.path.isdir(directory): os.system("mkdir "+directory)
for mA in mA_list:
  for mphi in mphi_list:
    filename = "phi{}A{}To4Tau".format(int(mphi),int(mA))
    if not os.path.isdir(directory+"/"+filename): os.system("mkdir "+directory+"/"+filename)

    # Determine alignment scenario
    if mphi > 125: # Normal alignment
      mh = 125
      mH = mphi
      sinbma = 1
      phi = "h2"
    else: # Inverted alignment
      mh = mphi
      mH = 125
      sinbma = 0
      phi = "h1"
    mixh   = (m.pi/2) - m.asin(sinbma)
    mHc = mphi

    find_list = ["__mA__","__mH__","__mh__","__mixh__","__mHc__","__phi__","__mphi__"]
    replace_list = [mA,mH,mh,mixh,mHc,phi,mphi]
    WriteListToFile(AppendLinesInList(customize_card,find_list,replace_list),directory+"/"+filename+"/"+filename+"_customizecards.dat")
    WriteListToFile(AppendLinesInList(proc_card,find_list,replace_list),directory+"/"+filename+"/"+filename+"_proc_card.dat")
    WriteListToFile(AppendLinesInList(run_card,find_list,replace_list),directory+"/"+filename+"/"+filename+"_run_card.dat")

    if args.generate and not args.submit: os.system("./gridpack_generation.sh {} {}/{}".format(filename,directory,filename))
    if args.submit:
      WriteListToFile(AppendLinesInList(submit_file,["__filename__","__directory__"],[filename,directory]),filename+"_gridpack_submissions.sh")
      os.system("qsub -e {}_gridpack_submissions_error.txt -o {}_gridpack_submissions_output.txt -V -q hep.q -l h_rt=0:600:0 -cwd {}_gridpack_submissions.sh".format(filename,filename,filename))
