import subprocess
from datetime import datetime

# Get the current date and time
current_datetime = datetime.now().strftime("%Y%m%d%H%M%S")

# main_input_dir = "/media/jfklamka/SAMSUNG_FUW/ILD/data"
main_input_dir = "/media/jfklamka/SAMSUNG_FUW/ILD/data/vertex_tests"
# main_input_dir = "/home/jfklamka/ILC/TrackingAnalysis/data/vertex_tests"

# Construct the output filename with date and time
output_filename = f"Efficiency/fimp_E250_simulation1_l_overlay_{current_datetime}_NdfCut_d0z0Cuts_isoCuts_R160cmZ2m_cosOpen0995.txt"
# output_filename = f"Efficiency/fimp_E250_simulation1_l_overlay_{current_datetime}_testRelMomentumCut.txt"

masses = [(109,110), (55,110), (12,110), (59,60), (30,60), (12,60)]
# masses = [ (59,60),(30,60)]

# Path to your main program
analysis = "fimp_analysis.py"

for m in masses:

    t = 3000
    # if m[1] == 60 or m == (109, 110) or m == (12, 110):
    if m[1] == 60:
        t = 1000

    # output_filename = f"outHistFiles/kink_tests/fimp_simulation1_l_Ms{m[0]}_M{m[1]}_overlay_NdfCut_d0z0Cuts_isoCuts_R160cmZ2m_cosOpen0995.root"
    output_filename = f"outHistFiles/kink_tests/fimp_simulation1_l_Ms{m[0]}_M{m[1]}_overlay_NdfCut.root"

    # input_filename = f"{main_input_dir}/fimp_simulation1_l_Ms{m[0]}*ev_M{m[1]}gev_t*_overlay_*.slcio"
    # input_filename = f"{main_input_dir}/fimp_simulation1_l_Ms{m[0]}*ev_M{m[1]}gev_t{t}_overlay_testRelMomentumCut.slcio"
    input_filename = f"{main_input_dir}/fimp_simulation1_l_Ms{m[0]}*ev_M{m[1]}gev_t{t}_overlay_R60mm_defaultHadronCuts_TPCcut30.slcio"

    # Construct the command to run the main program with the current configuration
    # command = ["python", analysis, input_filename, output_filename]
    command = "python " + analysis + " " + input_filename + " " + output_filename

    # Run the main program with subprocess
    subprocess.run(command, shell=True, check=True)