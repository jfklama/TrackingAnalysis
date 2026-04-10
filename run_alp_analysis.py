import subprocess
from datetime import datetime

# Get the current date and time
current_datetime = datetime.now().strftime("%Y%m%d%H%M%S")

main_input_dir = "/media/jfklamka/SAMSUNG_FUW/ILD/data/vertex_tests"
# main_input_dir = "/home/jfklamka/ILC/TrackingAnalysis/data/vertex_tests"

# Construct the output filename with date and time
output_filename = f"Efficiency/alp_test_simulation5_6_noFinalCuts_{current_datetime}.txt"

# masses = [300, 1000, 3000, 10000]
masses = [3000, 10000]

# Path to main program
analysis = "idm_analysis.py"

for m in masses:

    output_filename = f"outHistFiles/vertex_selection/alp_test_simulation5_6_M{m}_overlay_noCuts.root"
    # output_filename = f"outHistFiles/vertex_selection/alp_test_simulation5_6_M{m}_overlay_allCuts_refDist50_test.root"

    # input_filename = f"{main_input_dir}/alp_test_simulation5_6_M{m}_overlay_allCuts_refDist50.slcio"
    # input_filename = f"{main_input_dir}/alp_test_simulation5_6_M{m}_overlay_allCuts_allR.slcio"
    input_filename = f"{main_input_dir}/alp_test_simulation5_6_M{m}_overlay_noCuts.slcio"

    # Construct the command to run the main program with the current configuration
    command = ["python", analysis, input_filename, output_filename]

    # Run the main program with subprocess
    subprocess.run(command)