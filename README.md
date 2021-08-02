# TrackingAnalysis
 Analysis of the ILD tracking performance for LLP detection
 
 ## Examples
 
 ### SIMULATION

simulate events with particle gun:

```shell
ddsim --steeringFile ddsim_steer.py --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --enableGun --gun.particle lambda --gun.energy 20*GeV --gun.distribution uniform --gun.isotrop True -N 1000 --outputFile <outputFileName>.slcio
```

simulate from input samples (with particle table):

```shell
ddsim \
    --inputFiles data/idm_test_generation_dM_10.stdhep \
    -N 10000 \
    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
    --physics.pdgfile particle_idm_dM_10.tbl \
    --outputFile data/idm_test_simulation_dM_10.slcio
```

### RECONSTRUCTION

```shell
Marlin MarlinStdReco.xml --constant.lcgeo_DIR=$lcgeo_DIR --constant.DetectorModel=ILD_l5_o1_v02 --constant.OutputBaseName=<outputBaseName> --global.LCIOInputFiles=<inputFileName>.slcio
```
