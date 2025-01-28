
# 5 GeV
#ddsim --steeringFile ddsim_steer.py --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --enableGun --gun.particle lambda --gun.energy 5*GeV --gun.distribution uniform --gun.isotrop True -N 5000 --outputFile data/sm_lambdas/lambda_5gev_5000evt_gun.slcio

Marlin MarlinStdReco.xml --constant.lcgeo_DIR=$lcgeo_DIR --constant.DetectorModel=ILD_l5_o1_v02 --constant.OutputBaseName=data/sm_lambdas/lambda_5gev_5000evt --global.LCIOInputFiles=data/sm_lambdas/lambda_5gev_5000evt_gun.slcio

rm -r ./data/sm_lambdas/*REC.slcio

# 10 GeV
#ddsim --steeringFile ddsim_steer.py --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --enableGun --gun.particle lambda --gun.energy 10*GeV --gun.distribution uniform --gun.isotrop True -N 5000 --outputFile data/sm_lambdas/lambda_10gev_5000evt_gun.slcio

Marlin MarlinStdReco.xml --constant.lcgeo_DIR=$lcgeo_DIR --constant.DetectorModel=ILD_l5_o1_v02 --constant.OutputBaseName=data/sm_lambdas/lambda_10gev_5000evt --global.LCIOInputFiles=data/sm_lambdas/lambda_10gev_5000evt_gun.slcio

rm -r ./data/sm_lambdas/*REC.slcio

# 20 GeV
#ddsim --steeringFile ddsim_steer.py --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --enableGun --gun.particle lambda --gun.energy 20*GeV --gun.distribution uniform --gun.isotrop True -N 5000 --outputFile data/sm_lambdas/lambda_20gev_5000evt_gun.slcio

Marlin MarlinStdReco.xml --constant.lcgeo_DIR=$lcgeo_DIR --constant.DetectorModel=ILD_l5_o1_v02 --constant.OutputBaseName=data/sm_lambdas/lambda_20gev_5000evt --global.LCIOInputFiles=data/sm_lambdas/lambda_20gev_5000evt_gun.slcio

rm -r ./data/sm_lambdas/*REC.slcio

# 50 GeV
#ddsim --steeringFile ddsim_steer.py --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml --enableGun --gun.particle lambda --gun.energy 50*GeV --gun.distribution uniform --gun.isotrop True -N 5000 --outputFile data/sm_lambdas/lambda_50gev_5000evt_gun.slcio

Marlin MarlinStdReco.xml --constant.lcgeo_DIR=$lcgeo_DIR --constant.DetectorModel=ILD_l5_o1_v02 --constant.OutputBaseName=data/sm_lambdas/lambda_50gev_5000evt --global.LCIOInputFiles=data/sm_lambdas/lambda_50gev_5000evt_gun.slcio

rm -r ./data/sm_lambdas/*REC.slcio
