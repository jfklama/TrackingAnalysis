declare -a particles=(6) # 4=pi, 5=e, 6=mu, 7=tau

overlay_ee=$(ls /home/jfklamka/ILC/overlay/hep.fuw.edu.pl/u/zarnecki/ilc/overlay/00010233/000/*.slcio)
overlay_ww=$(ls /home/jfklamka/ILC/overlay/BgOverlayWW/*.slcio)
overlay_wb=$(ls /home/jfklamka/ILC/overlay/BgOverlayWB/*.slcio)
overlay_bw=$(ls /home/jfklamka/ILC/overlay/BgOverlayBW/*.slcio)
overlay_bb=$(ls /home/jfklamka/ILC/overlay/BgOverlayBB/*.slcio)


for p in ${particles[@]}
do

  #echo $overlay_files

  #ddsim --steeringFile ddsim_steer.py \
  #    --inputFiles data/idm_test_generation5_${p}.slcio \
  #    -N 1000 \
  #    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
  #    --outputFile data/idm_test_simulation5_${p}_dM50_SiliconTracker_1000ev.slcio

  # ddsim --steeringFile ddsim_steer.py \
  #    --inputFiles data/idm_test_generation5_${p}_dM20.slcio \
  #    -N 100000 \
  #    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
  #    --outputFile data/idm_test_simulation5_${p}_dM20.slcio \
  #    --physics.pdgfile ParticleLists/particle_idm_dM_20.tbl

  # ddsim --steeringFile ddsim_steer.py \
  #    --inputFiles data/idm_test_generation5_${p}_dM100.slcio \
  #    -N 100000 \
  #    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
  #    --outputFile data/idm_test_simulation5_${p}_dM100.slcio \
  #    --physics.pdgfile ParticleLists/particle_idm_dM_100.tbl


  Marlin MarlinStdReco.xml \
         --constant.lcgeo_DIR=$lcgeo_DIR \
         --constant.DetectorModel=ILD_l5_o1_v02 \
         --constant.OutputBaseName=data/idm_test_simulation5_${p}_dM10_noOverlay_100ev \
         --constant.CMSEnergy=500 \
         --global.LCIOInputFiles=data/idm_test_simulation5_${p}_dM10.slcio \
         --global.MaxRecordNumber=100 \
         --MyClupatraProcessor.CreateDebugCollections=true \
        #  --global.Verbosity=DEBUG \
         # --MyLLPFinder.CutOnHelixDistance=25 \
         # --MyLLPFinder.CutOnTrackCurvatureRatio=0.95 \
         # --MyLLPFinder.CutOnTrackOpeningAngle=1.0 \
         # --MyLLPFinder.CutOnVertexDistance=0 \
         # --constant.RunOverlay=true \
         # --PairBgOverlay.InputFileNames=$overlay_ee \
         # --BgOverlayWW.InputFileNames=$overlay_ww \
         # --BgOverlayWB.InputFileNames=$overlay_wb \
         # --BgOverlayBW.InputFileNames=$overlay_bw \
         # --BgOverlayBB.InputFileNames=$overlay_bb \
         # --constant.ExpectedBgWW=0.211 \
         # --constant.ExpectedBgWB=0.24605 \
         # --constant.ExpectedBgBW=0.243873 \
         # --constant.ExpectedBgBB=0.35063 \
         #--constant.RunOverlay=true

# Marlin MarlinStdReco.xml \
#          --constant.lcgeo_DIR=$lcgeo_DIR \
#          --constant.DetectorModel=ILD_l5_o1_v02 \
#          --constant.OutputBaseName=data/idm_test_simulation5_${p}_dM20_noOverlay_test \
#          --global.LCIOInputFiles=/media/jfklamka/Filip_USB3/ILC/data/idm_test_simulation5_${p}_dM20.slcio \
         # --MyLLPFinder.CutOnHelixDistance=25 \
         # --MyLLPFinder.CutOnTrackCurvatureRatio=0.95 \
         # --MyLLPFinder.CutOnTrackOpeningAngle=1.0 \
         # --MyLLPFinder.CutOnVertexDistance=0 \
         # --constant.RunOverlay=true \
         # --PairBgOverlay.InputFileNames=$overlay_ee \
         # --BgOverlayWW.InputFileNames=$overlay_ww \
         # --BgOverlayWB.InputFileNames=$overlay_wb \
         # --BgOverlayBW.InputFileNames=$overlay_bw \
         # --BgOverlayBB.InputFileNames=$overlay_bb \
         # --constant.ExpectedBgWW=0.211 \
         # --constant.ExpectedBgWB=0.24605 \
         # --constant.ExpectedBgBW=0.243873 \
         # --constant.ExpectedBgBB=0.35063 \
         #--global.MaxRecordNumber=1000 \
         #--global.Verbosity=DEBUG
         #--constant.RunOverlay=true

# Marlin MarlinStdReco.xml \
#          --constant.lcgeo_DIR=$lcgeo_DIR \
#          --constant.DetectorModel=ILD_l5_o1_v02 \
#          --constant.OutputBaseName=data/idm_test_simulation5_${p}_dM10_noOverlay_test \
#          --global.LCIOInputFiles=data/idm_test_simulation5_${p}_dM10.slcio \
         # --MyLLPFinder.CutOnHelixDistance=25 \
         # --MyLLPFinder.CutOnTrackCurvatureRatio=0.95 \
         # --MyLLPFinder.CutOnTrackOpeningAngle=1.0 \
         # --MyLLPFinder.CutOnVertexDistance=0 \
         # --constant.RunOverlay=true \
         # --PairBgOverlay.InputFileNames=$overlay_ee \
         # --BgOverlayWW.InputFileNames=$overlay_ww \
         # --BgOverlayWB.InputFileNames=$overlay_wb \
         # --BgOverlayBW.InputFileNames=$overlay_bw \
         # --BgOverlayBB.InputFileNames=$overlay_bb \
         # --constant.ExpectedBgWW=0.211 \
         # --constant.ExpectedBgWB=0.24605 \
         # --constant.ExpectedBgBW=0.243873 \
         # --constant.ExpectedBgBB=0.35063 \
         #--global.MaxRecordNumber=1000 \
         #--global.Verbosity=DEBUG
         #--constant.RunOverlay=true

# Marlin MarlinStdReco.xml \
#          --constant.lcgeo_DIR=$lcgeo_DIR \
#          --constant.DetectorModel=ILD_l5_o1_v02 \
#          --constant.OutputBaseName=data/idm_test_simulation5_${p}_dM30_noOverlay_test \
#          --global.LCIOInputFiles=/media/jfklamka/Filip_USB3/ILC/data/idm_test_simulation5_${p}_dM30.slcio \
#          # --MyLLPFinder.CutOnHelixDistance=25 \
#          # --MyLLPFinder.CutOnTrackCurvatureRatio=0.95 \
#          # --MyLLPFinder.CutOnTrackOpeningAngle=1.0 \
#          # --MyLLPFinder.CutOnVertexDistance=0 \
#          # --constant.RunOverlay=true \
#          # --PairBgOverlay.InputFileNames=$overlay_ee \
#          # --BgOverlayWW.InputFileNames=$overlay_ww \
#          # --BgOverlayWB.InputFileNames=$overlay_wb \
#          # --BgOverlayBW.InputFileNames=$overlay_bw \
#          # --BgOverlayBB.InputFileNames=$overlay_bb \
#          # --constant.ExpectedBgWW=0.211 \
#          # --constant.ExpectedBgWB=0.24605 \
#          # --constant.ExpectedBgBW=0.243873 \
#          # --constant.ExpectedBgBB=0.35063 \
#          # --global.MaxRecordNumber=1000 \
#          #--global.Verbosity=DEBUG
#          #--constant.RunOverlay=true

# Marlin MarlinStdReco.xml \
#          --constant.lcgeo_DIR=$lcgeo_DIR \
#          --constant.DetectorModel=ILD_l5_o1_v02 \
#          --constant.OutputBaseName=data/idm_test_simulation5_${p}_dM100_noOverlay_1000ev \
#          --global.LCIOInputFiles=data/idm_test_simulation5_${p}_dM100.slcio \
#          --global.MaxRecordNumber=500 \
#          # --MyLLPFinder.CutOnHelixDistance=25 \
#          # --MyLLPFinder.CutOnTrackCurvatureRatio=0.95 \
#          # --MyLLPFinder.CutOnTrackOpeningAngle=1.0 \
#          # --MyLLPFinder.CutOnVertexDistance=0 \
#          # --constant.RunOverlay=true \
#          # --PairBgOverlay.InputFileNames=$overlay_ee \
#          # --BgOverlayWW.InputFileNames=$overlay_ww \
#          # --BgOverlayWB.InputFileNames=$overlay_wb \
#          # --BgOverlayBW.InputFileNames=$overlay_bw \
#          # --BgOverlayBB.InputFileNames=$overlay_bb \
#          # --constant.ExpectedBgWW=0.211 \
#          # --constant.ExpectedBgWB=0.24605 \
#          # --constant.ExpectedBgBW=0.243873 \
#          # --constant.ExpectedBgBB=0.35063 \
#          #--global.Verbosity=DEBUG
#          #--constant.RunOverlay=true


done
