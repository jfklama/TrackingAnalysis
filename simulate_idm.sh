
# dM = 1 GeV
ddsim \
    --inputFiles data/idm_test_generation_dM_10.stdhep \
    -N 10000 \
    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
    --physics.pdgfile ParticleLists/particle_idm_dM_10.tbl \
    --outputFile data/idm_test_simulation_dM_10.slcio

Marlin MarlinStdReco.xml \
       --constant.lcgeo_DIR=$lcgeo_DIR \
       --constant.DetectorModel=ILD_l5_o1_v02 \
       --constant.OutputBaseName=data/idm_test_simulation_dM_10 \
       --global.LCIOInputFiles=data/idm_test_simulation_dM_10.slcio

# dM = 5 GeV
ddsim \
    --inputFiles data/idm_test_generation_dM_50.stdhep \
    -N 10000 \
    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
    --physics.pdgfile ParticleLists/particle_idm_dM_50.tbl \
    --outputFile data/idm_test_simulation_dM_50.slcio

Marlin MarlinStdReco.xml \
       --constant.lcgeo_DIR=$lcgeo_DIR \
       --constant.DetectorModel=ILD_l5_o1_v02 \
       --constant.OutputBaseName=data/idm_test_simulation_dM_50 \
       --global.LCIOInputFiles=data/idm_test_simulation_dM_50.slcio
