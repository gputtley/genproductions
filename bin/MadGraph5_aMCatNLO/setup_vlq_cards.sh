betaRd33_0_betaL32=-0.15
betaRd33_0_betaL32_1up=-0.02
betaRd33_0_betaL32_1down=-0.26
betaRd33_0_betaL32_2up=0.11
betaRd33_0_betaL32_2down=-0.37


betaRd33_0_betaL23=0.19
betaRd33_0_betaL23_1up=0.25
betaRd33_0_betaL23_1down=0.10
betaRd33_0_betaL23_2up=0.0.31
betaRd33_0_betaL23_2down=0.0.01


betaRd33_minus1_betaL32=-0.21
betaRd33_minus1_betaL32_1up=-0.16
betaRd33_minus1_betaL32_1down=-0.26
betaRd33_minus1_betaL32_2up=-0.11
betaRd33_minus1_betaL32_2down=-0.31

betaRd33_minus1_betaL23=0.21
betaRd33_minus1_betaL23_1up=0.26
betaRd33_minus1_betaL23_1down=0.12
betaRd33_minus1_betaL23_2up=0.31
betaRd33_minus1_betaL23_2down=0.02



mkdir cards/vlq

declare -A mU_gU

mU_gU[0,0]="2"
mU_gU[0,1]="1"

mU_gU[1,0]="2"
mU_gU[1,1]="2"

mU_gU[2,0]="2"
mU_gU[2,1]="3"

mU_gU[3,0]="3"
mU_gU[3,1]="1"

mU_gU[4,0]="3"
mU_gU[4,1]="2"

mU_gU[5,0]="3"
mU_gU[5,1]="3"

mU_gU[6,0]="4"
mU_gU[6,1]="1"

mU_gU[7,0]="4"
mU_gU[7,1]="2"

mU_gU[8,0]="4"
mU_gU[8,1]="3"


for ((j=0;j<=8;j++)) 
do
  for i in betaRd33_0 betaRd33_minus1
  do

    betaL23_string="${i}_betaL23"
    betaL23_1up_string="${i}_betaL23_1up"
    betaL23_1down_string="${i}_betaL23_1down"
    betaL23_2up_string="${i}_betaL23_2up"
    betaL23_2down_string="${i}_betaL23_2down"

    betaL32_string="${i}_betaL32"
    betaL32_1up_string="${i}_betaL32_1up"
    betaL32_1down_string="${i}_betaL32_1down"
    betaL32_2up_string="${i}_betaL32_2up"
    betaL32_2down_string="${i}_betaL32_2down"

    mkdir "cards/vlq/${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_}"
    filename="cards/vlq/${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_}/${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_}"
    
    # Set up customizecard
    echo "set param_card mass 9000007 ${mU_gU[$j,0]}" > "${filename}_customizecards.dat"
    echo "set param_card nplqcoup 1 ${mU_gU[$j,1]}" >> "${filename}_customizecards.dat"
    echo "set param_card nplqcoup 2 1.000000e+00" >> "${filename}_customizecards.dat"

    if [ $i == "betaRd33_0" ]
    then
      echo "set param_card nplqcoup 3 0.000000e+00" >> "${filename}_customizecards.dat"
      echo "set param_card nplqcoup 4 ${betaRd33_0_betaL23}" >> "${filename}_customizecards.dat"
      echo "set param_card nplqcoup 5 ${betaRd33_0_betaL32}" >> "${filename}_customizecards.dat"
    elif [ $i == "betaRd33_minus1" ]
    then
      echo "set param_card nplqcoup 3 -1.000000e+00" >> "${filename}_customizecards.dat"
      echo "set param_card nplqcoup 4 ${betaRd33_minus1_betaL23}" >> "${filename}_customizecards.dat"
      echo "set param_card nplqcoup 5 ${betaRd33_minus1_betaL32}" >> "${filename}_customizecards.dat"
    fi    

    # Set up proc_card
    echo "set default_unset_couplings 99" > "${filename}_proc_card.dat"
    echo "set group_subprocesses Auto" >> "${filename}_proc_card.dat"
    echo "set ignore_six_quark_processes False" >> "${filename}_proc_card.dat"
    echo "set loop_optimized_output True" >> "${filename}_proc_card.dat"
    echo "set loop_color_flows False" >> "${filename}_proc_card.dat"
    echo "set gauge unitary" >> "${filename}_proc_card.dat"
    echo "set complex_mass_scheme False" >> "${filename}_proc_card.dat"
    echo "set max_npoint_for_channel 0" >> "${filename}_proc_card.dat"
    echo "import model vector_LQ_UFO" >> "${filename}_proc_card.dat"
    echo "define p = 21 2 4 1 3 -2 -4 -1 -3 5 -5 # pass to 5 flavors" >> "${filename}_proc_card.dat"
    echo "define j = p" >> "${filename}_proc_card.dat"
    echo "generate p p > ta+ ta- / zp gp z a" >> "${filename}_proc_card.dat"
    echo "output ${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_} -nojpeg" >> "${filename}_proc_card.dat"


    # Set up run_card
    echo "#*********************************************************************" > "${filename}_run_card.dat"
    echo "#                       MadGraph5_aMC@NLO                            *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#                     run_card.dat MadEvent                          *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#  This file is used to set the parameters of the run.               *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#  Some notation/conventions:                                        *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#   Lines starting with a '# ' are info or comments                  *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#   mind the format:   value    = variable     ! comment             *" >> "${filename}_run_card.dat"
    echo "#                                                                    *" >> "${filename}_run_card.dat"
    echo "#   To display more options, you can type the command:               *" >> "${filename}_run_card.dat"
    echo "#      update full_run_card                                          *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "#*******************" >> "${filename}_run_card.dat"                                                 
    echo "# Running parameters" >> "${filename}_run_card.dat"
    echo "#*******************" >> "${filename}_run_card.dat"                                                 
    echo "#" >> "${filename}_run_card.dat"                                                                    
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Tag name for the run (one word)                                    *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "  tag_1     = run_tag ! name of the run" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Number of events and rnd seed                                      *" >> "${filename}_run_card.dat"
    echo "# Warning: Do not generate more than 1M events in a single run       *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "  10000 = nevents ! Number of unweighted events requested" >> "${filename}_run_card.dat"
    echo "  0   = iseed   ! rnd seed (0=assigned automatically=default))" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Collider type and energy                                           *" >> "${filename}_run_card.dat"
    echo "# lpp: 0=No PDF, 1=proton, -1=antiproton, 2=photon from proton,      *" >> "${filename}_run_card.dat"
    echo "#                                         3=photon from electron     *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "     1        = lpp1    ! beam 1 type" >> "${filename}_run_card.dat"
    echo "     1        = lpp2    ! beam 2 type" >> "${filename}_run_card.dat"
    echo "     6500.0     = ebeam1  ! beam 1 total energy in GeV" >> "${filename}_run_card.dat"
    echo "     6500.0     = ebeam2  ! beam 2 total energy in GeV" >> "${filename}_run_card.dat"
    echo "# To see polarised beam options: type "update beam_pol"" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# PDF CHOICE: this automatically fixes also alpha_s and its evol.    *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "     lhapdf = pdlabel ! PDF set" >> "${filename}_run_card.dat"
    echo "     \$DEFAULT_PDF_SETS = lhaid ! If pdlabel=lhapdf, this is the lhapdf number. Only" >> "${filename}_run_card.dat"
    echo "                  ! numbers for central PDF sets are allowed. Can be a list;" >> "${filename}_run_card.dat"
    echo "                  ! PDF sets beyond the first are included via reweighting." >> "${filename}_run_card.dat"
    echo "# To see heavy ion options: type \"update ion_pdf\" " >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Renormalization and factorization scales                           *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " False = fixed_ren_scale  ! if .true. use fixed ren scale" >> "${filename}_run_card.dat"
    echo " False        = fixed_fac_scale  ! if .true. use fixed fac scale" >> "${filename}_run_card.dat"
    echo " 91.188  = scale            ! fixed ren scale" >> "${filename}_run_card.dat"
    echo " 91.188  = dsqrt_q2fact1    ! fixed fact scale for pdf1" >> "${filename}_run_card.dat"
    echo " 91.188  = dsqrt_q2fact2    ! fixed fact scale for pdf2" >> "${filename}_run_card.dat"
    echo " -1 = dynamical_scale_choice ! Choose one of the preselected dynamical choices" >> "${filename}_run_card.dat"
    echo " 1.0  = scalefact        ! scale factor for event-by-event scales" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Type and output format" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "  False     = gridpack  !True = setting up the grid pack" >> "${filename}_run_card.dat"
    echo "  -1.0 = time_of_flight ! threshold (in mm) below which the invariant livetime is not written (-1 means not written)" >> "${filename}_run_card.dat"
    echo "  3.0 = lhe_version       ! Change the way clustering information pass to shower." >> "${filename}_run_card.dat"
    echo "  True = clusinfo         ! include clustering tag in output" >> "${filename}_run_card.dat"
    echo "  average =  event_norm       ! average/sum. Normalization of the weight in the LHEF" >> "${filename}_run_card.dat"
    echo "" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Matching parameter (MLM only)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0 = ickkw            ! 0 no matching, 1 MLM" >> "${filename}_run_card.dat"
    echo " 1.0 = alpsfact         ! scale factor for QCD emission vx" >> "${filename}_run_card.dat"
    echo " False = chcluster        ! cluster only according to channel diag" >> "${filename}_run_card.dat"
    echo " 5 = asrwgtflavor     ! highest quark flavor for a_s reweight" >> "${filename}_run_card.dat"
    echo " True  = auto_ptj_mjj  ! Automatic setting of ptj and mjj if xqcut >0" >> "${filename}_run_card.dat"
    echo "                                   ! (turn off for VBF and single top processes)" >> "${filename}_run_card.dat"
    echo " 0.0   = xqcut   ! minimum kt jet measure between partons" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# handling of the helicities:" >> "${filename}_run_card.dat"
    echo "#  0: sum over all helicities" >> "${filename}_run_card.dat"
    echo "#  1: importance sampling over helicities" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "   0  = nhel          ! using helicities importance sampling or not." >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Generation bias, check the wiki page below for more information:   *" >> "${filename}_run_card.dat"
    echo "#  'cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias' *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " None = bias_module  ! Bias type of bias, [None, ptj_bias, -custom_folder-]" >> "${filename}_run_card.dat"
    echo " {} = bias_parameters ! Specifies the parameters of the module." >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "#*******************************" >> "${filename}_run_card.dat"                                                 
    echo "# Parton level cuts definition *" >> "${filename}_run_card.dat"
    echo "#*******************************" >> "${filename}_run_card.dat"                                     
    echo "#" >> "${filename}_run_card.dat"                                                                    
    echo "#" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# BW cutoff (M+/-bwcutoff*Gamma) ! Define on/off-shell for \"$\" and decay" >> "${filename}_run_card.dat"  
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "  15.0  = bwcutoff      ! (M+/-bwcutoff*Gamma)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Apply pt/E/eta/dr/mij/kt_durham cuts on decay products or not" >> "${filename}_run_card.dat"
    echo "# (note that etmiss/ptll/ptheavy/ht/sorted cuts always apply)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "   False  = cut_decays    ! Cut decay products" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Standard Cuts                                                      *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Minimum and maximum pt's (for max, -1 means no cut)                *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 20.0  = ptj       ! minimum pt for the jets" >> "${filename}_run_card.dat"
    echo " 0.0  = ptb       ! minimum pt for the b" >> "${filename}_run_card.dat"
    echo " 10.0  = pta       ! minimum pt for the photons" >> "${filename}_run_card.dat"
    echo " 10.0  = ptl       ! minimum pt for the charged leptons" >> "${filename}_run_card.dat"
    echo " 0.0  = misset    ! minimum missing Et (sum of neutrino's momenta)" >> "${filename}_run_card.dat"
    echo " -1.0  = ptjmax    ! maximum pt for the jets" >> "${filename}_run_card.dat"
    echo " -1.0  = ptbmax    ! maximum pt for the b" >> "${filename}_run_card.dat"
    echo " -1.0  = ptamax    ! maximum pt for the photons" >> "${filename}_run_card.dat"
    echo " -1.0  = ptlmax    ! maximum pt for the charged leptons" >> "${filename}_run_card.dat"
    echo " -1.0  = missetmax ! maximum missing Et (sum of neutrino's momenta)" >> "${filename}_run_card.dat"
    echo " {} = pt_min_pdg ! pt cut for other particles (use pdg code). Applied on particle and anti-particle" >> "${filename}_run_card.dat"
    echo " {}     = pt_max_pdg ! pt cut for other particles (syntax e.g. {6: 100, 25: 50})" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "# For display option for energy cut in the partonic center of mass frame type 'update ecut'" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Maximum and minimum absolute rapidity (for max, -1 means no cut)   *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "  5.0 = etaj    ! max rap for the jets" >> "${filename}_run_card.dat"
    echo "  -1.0  = etab    ! max rap for the b" >> "${filename}_run_card.dat"
    echo " 2.5  = etaa    ! max rap for the photons" >> "${filename}_run_card.dat"
    echo " 2.5  = etal    ! max rap for the charged leptons" >> "${filename}_run_card.dat"
    echo " 0.0  = etajmin ! min rap for the jets" >> "${filename}_run_card.dat"
    echo " 0.0  = etabmin ! min rap for the b" >> "${filename}_run_card.dat"
    echo " 0.0  = etaamin ! min rap for the photons" >> "${filename}_run_card.dat"
    echo " 0.0  = etalmin ! main rap for the charged leptons" >> "${filename}_run_card.dat"
    echo " {} = eta_min_pdg ! rap cut for other particles (use pdg code). Applied on particle and anti-particle" >> "${filename}_run_card.dat"
    echo " {} = eta_max_pdg ! rap cut for other particles (syntax e.g. {6: 2.5, 23: 5})" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Minimum and maximum DeltaR distance                                *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.4 = drjj    ! min distance between jets" >> "${filename}_run_card.dat"
    echo " 0.0   = drbb    ! min distance between b's" >> "${filename}_run_card.dat"
    echo " 0.4 = drll    ! min distance between leptons" >> "${filename}_run_card.dat"
    echo " 0.4 = draa    ! min distance between gammas" >> "${filename}_run_card.dat"
    echo " 0.0   = drbj    ! min distance between b and jet" >> "${filename}_run_card.dat"
    echo " 0.4 = draj    ! min distance between gamma and jet" >> "${filename}_run_card.dat"
    echo " 0.4 = drjl    ! min distance between jet and lepton" >> "${filename}_run_card.dat"
    echo " 0.0   = drab    ! min distance between gamma and b" >> "${filename}_run_card.dat"
    echo " 0.0   = drbl    ! min distance between b and lepton" >> "${filename}_run_card.dat"
    echo " 0.4 = dral    ! min distance between gamma and lepton" >> "${filename}_run_card.dat"
    echo " -1.0  = drjjmax ! max distance between jets" >> "${filename}_run_card.dat"
    echo " -1.0  = drbbmax ! max distance between b's" >> "${filename}_run_card.dat"
    echo " -1.0  = drllmax ! max distance between leptons" >> "${filename}_run_card.dat"
    echo " -1.0  = draamax ! max distance between gammas" >> "${filename}_run_card.dat"
    echo " -1.0  = drbjmax ! max distance between b and jet" >> "${filename}_run_card.dat"
    echo " -1.0  = drajmax ! max distance between gamma and jet" >> "${filename}_run_card.dat"
    echo " -1.0  = drjlmax ! max distance between jet and lepton" >> "${filename}_run_card.dat"
    echo " -1.0  = drabmax ! max distance between gamma and b" >> "${filename}_run_card.dat"
    echo " -1.0  = drblmax ! max distance between b and lepton" >> "${filename}_run_card.dat"
    echo " -1.0  = dralmax ! maxdistance between gamma and lepton" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Minimum and maximum invariant mass for pairs                       *" >> "${filename}_run_card.dat"
    echo "# WARNING: for four lepton final state mmll cut require to have      *" >> "${filename}_run_card.dat"
    echo "#          different lepton masses for each flavor!                  *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = mmjj    ! min invariant mass of a jet pair" >> "${filename}_run_card.dat"
    echo " 0.0   = mmbb    ! min invariant mass of a b pair" >> "${filename}_run_card.dat"
    echo " 0.0   = mmaa    ! min invariant mass of gamma gamma pair" >> "${filename}_run_card.dat"
    echo " 0.0   = mmll    ! min invariant mass of l+l- (same flavour) lepton pair" >> "${filename}_run_card.dat"
    echo " -1.0  = mmjjmax ! max invariant mass of a jet pair" >> "${filename}_run_card.dat"
    echo " -1.0  = mmbbmax ! max invariant mass of a b pair" >> "${filename}_run_card.dat"
    echo " -1.0  = mmaamax ! max invariant mass of gamma gamma pair" >> "${filename}_run_card.dat"
    echo " -1.0  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair" >> "${filename}_run_card.dat"
    echo " {} = mxx_min_pdg ! min invariant mass of a pair of particles X/X~ (e.g. {6:250})" >> "${filename}_run_card.dat"
    echo " {'default': False} = mxx_only_part_antipart ! if True the invariant mass is applied only" >> "${filename}_run_card.dat"
    echo "                       ! to pairs of particle/antiparticle and not to pairs of the same pdg codes." >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Minimum and maximum invariant mass for all letpons                 *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = mmnl    ! min invariant mass for all letpons (l+- and vl)" >> "${filename}_run_card.dat"
    echo " -1.0  = mmnlmax ! max invariant mass for all letpons (l+- and vl)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Minimum and maximum pt for 4-momenta sum of leptons                *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = ptllmin  ! Minimum pt for 4-momenta sum of leptons(l and vl)" >> "${filename}_run_card.dat"
    echo " -1.0  = ptllmax  ! Maximum pt for 4-momenta sum of leptons(l and vl)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Inclusive cuts                                                     *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0  = ptheavy   ! minimum pt for at least one heavy final state" >> "${filename}_run_card.dat"
    echo " 0.0  = xptj ! minimum pt for at least one jet" >> "${filename}_run_card.dat"
    echo " 0.0  = xptb ! minimum pt for at least one b" >> "${filename}_run_card.dat"
    echo " 0.0  = xpta ! minimum pt for at least one photon" >> "${filename}_run_card.dat"
    echo " 0.0  = xptl ! minimum pt for at least one charged lepton" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Control the pt's of the jets sorted by pt                          *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = ptj1min ! minimum pt for the leading jet in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptj2min ! minimum pt for the second jet in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptj3min ! minimum pt for the third jet in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptj4min ! minimum pt for the fourth jet in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptj1max ! maximum pt for the leading jet in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptj2max ! maximum pt for the second jet in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptj3max ! maximum pt for the third jet in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptj4max ! maximum pt for the fourth jet in pt" >> "${filename}_run_card.dat"
    echo " 0   = cutuse  ! reject event if fails any (0) / all (1) jet pt cuts" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Control the pt's of leptons sorted by pt                           *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = ptl1min ! minimum pt for the leading lepton in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptl2min ! minimum pt for the second lepton in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptl3min ! minimum pt for the third lepton in pt" >> "${filename}_run_card.dat"
    echo " 0.0   = ptl4min ! minimum pt for the fourth lepton in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptl1max ! maximum pt for the leading lepton in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptl2max ! maximum pt for the second lepton in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptl3max ! maximum pt for the third lepton in pt" >> "${filename}_run_card.dat"
    echo " -1.0  = ptl4max ! maximum pt for the fourth lepton in pt" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Control the Ht(k)=Sum of k leading jets                            *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = htjmin ! minimum jet HT=Sum(jet pt)" >> "${filename}_run_card.dat"
    echo " -1.0  = htjmax ! maximum jet HT=Sum(jet pt)" >> "${filename}_run_card.dat"
    echo " 0.0   = ihtmin  !inclusive Ht for all partons (including b)" >> "${filename}_run_card.dat"
    echo " -1.0  = ihtmax  !inclusive Ht for all partons (including b)" >> "${filename}_run_card.dat"
    echo " 0.0   = ht2min ! minimum Ht for the two leading jets" >> "${filename}_run_card.dat"
    echo " 0.0   = ht3min ! minimum Ht for the three leading jets" >> "${filename}_run_card.dat"
    echo " 0.0   = ht4min ! minimum Ht for the four leading jets" >> "${filename}_run_card.dat"
    echo " -1.0  = ht2max ! maximum Ht for the two leading jets" >> "${filename}_run_card.dat"
    echo " -1.0  = ht3max ! maximum Ht for the three leading jets" >> "${filename}_run_card.dat"
    echo " -1.0  = ht4max ! maximum Ht for the four leading jets" >> "${filename}_run_card.dat"
    echo "#***********************************************************************" >> "${filename}_run_card.dat"
    echo "# Photon-isolation cuts, according to hep-ph/9801442                   *" >> "${filename}_run_card.dat"
    echo "# When ptgmin=0, all the other parameters are ignored                  *" >> "${filename}_run_card.dat"
    echo "# When ptgmin>0, pta and draj are not going to be used                 *" >> "${filename}_run_card.dat"
    echo "#***********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0 = ptgmin ! Min photon transverse momentum" >> "${filename}_run_card.dat"
    echo " 0.4 = R0gamma ! Radius of isolation code" >> "${filename}_run_card.dat"
    echo " 1.0 = xn ! n parameter of eq.(3.4) in hep-ph/9801442" >> "${filename}_run_card.dat"
    echo " 1.0 = epsgamma ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442" >> "${filename}_run_card.dat"
    echo " True = isoEM ! isolate photons from EM energy (photons and leptons)" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# WBF cuts                                                           *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 0.0   = xetamin ! minimum rapidity for two jets in the WBF case" >> "${filename}_run_card.dat"
    echo " 0.0   = deltaeta ! minimum rapidity for two jets in the WBF case" >> "${filename}_run_card.dat"
    echo "#***********************************************************************" >> "${filename}_run_card.dat"
    echo "# Turn on either the ktdurham or ptlund cut to activate                *" >> "${filename}_run_card.dat"
    echo "# CKKW(L) merging with Pythia8 [arXiv:1410.3012, arXiv:1109.4829]      *" >> "${filename}_run_card.dat"
    echo "#***********************************************************************" >> "${filename}_run_card.dat"
    echo " -1.0  =  ktdurham" >> "${filename}_run_card.dat"
    echo " 0.4   =  dparameter" >> "${filename}_run_card.dat"
    echo " -1.0  =  ptlund" >> "${filename}_run_card.dat"
    echo " 1, 2, 3, 4, 5, 6, 21, 9000005, 9000007  =  pdgs_for_merging_cut ! PDGs for two cuts above" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# maximal pdg code for quark to be considered as a light jet         *" >> "${filename}_run_card.dat"
    echo "# (otherwise b cuts are applied)                                     *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo " 5 = maxjetflavor    ! Maximum jet pdg code" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "# Store info for systematics studies                                 *" >> "${filename}_run_card.dat"
    echo "# WARNING: Do not use for interference type of computation           *" >> "${filename}_run_card.dat"
    echo "#*********************************************************************" >> "${filename}_run_card.dat"
    echo "   True  = use_syst      ! Enable systematics studies" >> "${filename}_run_card.dat"
    echo "#" >> "${filename}_run_card.dat"
    echo "systematics = systematics_program ! none, systematics [python], SysCalc [depreceted, C++]" >> "${filename}_run_card.dat"
    echo "['--mur=0.5,1,2', '--muf=0.5,1,2', '--pdf=errorset'] = systematics_arguments ! see: https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Systematics#Systematicspythonmodule" >> "${filename}_run_card.dat"
    echo "# Syscalc is deprecated but to see the associate options type'update syscalc'" >> "${filename}_run_card.dat"

    # Set up reweight_card
    echo "change rwgt_dir ./rwgt" > "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=off_diag_0" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 0.000000e+00" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 0.000000e+00" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL23_1sigma_up" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_1up_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL23_1sigma_down" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_1down_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL32_1sigma_up" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_1up_string}" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL32_1sigma_down" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_1down_string}" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL23_2sigma_up" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_2up_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    echo "launch --rwgt_name=betaL23_2sigma_down" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 4 ${!betaL23_2down_string}" >> "${filename}_reweight_card.dat"
    echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    # won't include as XS doesn't depend on betaL32 that much
    #echo "launch --rwgt_name=betaL32_2sigma_up" >> "${filename}_reweight_card.dat"
    #echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
    #echo "  set nplqcoup 5 ${!betaL32_2up_string}" >> "${filename}_reweight_card.dat"
    #echo "launch --rwgt_name=betaL32_2sigma_down" >> "${filename}_reweight_card.dat"
    #echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
    #echo "  set nplqcoup 5 ${!betaL32_2down_string}" >> "${filename}_reweight_card.dat"
    if [ ${mU_gU[$j,1]} != 1 ]
    then
      echo "launch --rwgt_name=gU_1" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 1 1" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    fi
    if [ ${mU_gU[$j,1]} != 2 ]
    then
      echo "launch --rwgt_name=gU_2" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 1 2" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    fi
    if [ ${mU_gU[$j,1]} != 3 ]
    then
      echo "launch --rwgt_name=gU_3" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 1 3" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    fi
    if [ $i != "betaRd33_0" ]
    then
      echo "launch --rwgt_name=betaR33_0" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 1 ${mU_gU[$j,1]}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 3 0" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    fi
    if [ $i != "betaRd33_minus1" ]
    then
      echo "launch --rwgt_name=betaR33_minus1" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 1 ${mU_gU[$j,1]}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 3 -1" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 4 ${!betaL23_string}" >> "${filename}_reweight_card.dat"
      echo "  set nplqcoup 5 ${!betaL32_string}" >> "${filename}_reweight_card.dat"
    fi



   # produce gridpacks
   eval "./gridpack_generation.sh ${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_} cards/vlq/${i}_mU${mU_gU[$j,0]/'.'/_}_gU${mU_gU[$j,1]/'.'/_}"

  done
done











