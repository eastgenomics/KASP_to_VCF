import os
import sys
import glob
import subprocess


def main(genotyping_data_dir="/mnt/storage/samba/samba.ctrulab.uk/cytogenetics/kasp_genotyping/"):

    runfolders = next(os.walk(genotyping_data_dir))[1] # All subdirectories of dir

    print("Processing runfolders in %s" % genotyping_data_dir)
    for runfolder in sorted(runfolders):
        print()
        print("RUN:", runfolder)

        runfolder_path = os.path.join(os.path.abspath(genotyping_data_dir),runfolder)
        print("PATH:", runfolder_path)
        
        done_file_glob_pattern = os.path.join(runfolder_path, "*.done")
        runfolder_done = any(glob.glob(done_file_glob_pattern))

        if runfolder_done:
            print("Data already processed")

        else:
            kasp_data_glob_pattern = os.path.join(runfolder_path, "*.csv")
            kasp_data = glob.glob(kasp_data_glob_pattern)
            
            if not kasp_data:
                print("Data not found")

            else:
                # These are our unprocessed output files
                print("Data awaiting processing")
                assert len(kasp_data) == 1, "Error, multiple kasp data files detected"
                project_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
                script_path = os.path.join(project_dir, "bin/KASPcsv_2_GELvcf.py")
                command = "python3.4 {script_path} {kasp_data}".format(script_path=script_path, kasp_data=kasp_data[0])
                print(command)
                subprocess.call(command, shell=True)

def usage():
    print("Usage:\nconvert_all_KASP_output.py <genotyping_data_dir>")

if __name__ == "__main__":
    main()