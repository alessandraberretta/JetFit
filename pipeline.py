from download import download_data
from manipulate_txt import manipulate_txt
from cosm_calc import cosm_calc
import subprocess
from argparse import ArgumentParser
import os


def write_bash_script(file, redshift, dist_lum, args):

    with open('script.sh', 'w') as script:
        script.write("#!/usr/bin/env bash\n")
        script.write("source activate JetFit\n")
        if args.localrepo:
            if args.pathout:
                script.write(
                    f"python fitter.py --grb {file} --z {redshift} --dL {dist_lum} --localrepo {args.localrepo} --pathout {args.pathout}")
            else:
                script.write(
                    f"python fitter.py --grb {file} --z {redshift} --dL {dist_lum} --localrepo {args.localrepo}")
        else:
            if args.pathout:
                script.write(
                    f"python fitter.py --grb {file} --z {redshift} --dL {dist_lum} --pathout {args.pathout}")
            else:
                script.write(
                    f"python fitter.py --grb {file} --z {redshift} --dL {dist_lum}")
    subprocess.run("chmod +x script.sh", shell=True, check=True)


def write_bash_script_condor(file, redshift, dist_lum, job_dir, fitter_dir, args):

    with open(job_dir + '/script.sh', 'w') as script:
        script.write("#!/usr/bin/env bash\n")
        script.write(
            "source /cvmfs/fermi.local.repo/anaconda3/etc/profile.d/conda.sh\n")
        script.write("conda activate fermi\n")
        if args.localrepo:
            if args.pathout:
                script.write(
                    f"python {fitter_dir}/fitter.py --grb {file} --z {redshift} --dL {dist_lum} --localrepo {args.localrepo} --pathout {args.pathout}")
            else:
                script.write(
                    f"python {fitter_dir}/fitter.py --grb {file} --z {redshift} --dL {dist_lum} --localrepo {args.localrepo}")
        else:
            if args.pathout:
                script.write(
                    f"python {fitter_dir}/fitter.py --grb {file} --z {redshift} --dL {dist_lum} --pathout {args.pathout}")
            else:
                script.write(
                    f"python {fitter_dir}/fitter.py --grb {file} --z {redshift} --dL {dist_lum}")
    subprocess.run(f"chmod +x {job_dir}/script.sh", shell=True, check=True)


def write_file_sub(job_dir):

    with open(job_dir + '/sub', 'w') as script:
        script.write("universe   = vanilla\n")
        script.write(f"executable = {job_dir}/script.sh\n")
        script.write(f"log={job_dir}/output.log\n")
        script.write(f"output={job_dir}/output.out\n")
        script.write(f"error={job_dir}/output.error\n")
        script.write('+OWNER="AlessandraB"\n')
        script.write("queue")


def use_downloaded_data(pathdata):

    list_file = [pathdata + '/' +
                 file for file in os.listdir(pathdata) if file.endswith('.data')]
    list_redshift = [float(i[i.rfind('_')+1:i.rfind('.')]) for i in list_file]

    print(f"Has been find {len(list_file)} in data folder: {pathdata}")
    return list_file, list_redshift


def main(args=None):

    parser = ArgumentParser(description='Pipeline for GRB analysis')

    parser.add_argument("-grb", "--grb", type=str,
                        dest="grb_name", help="GRB name")
    parser.add_argument("-year", "--year", type=int,
                        dest="grb_year", help="GRB year")
    parser.add_argument("-condor", "--condor",
                        dest="condor", action='store_true', help="condor")
    parser.add_argument("-pathdir", "--pathdir", type=str,
                        dest="pathdir", help="Condor job directory")
    parser.add_argument("-pathdata", "--pathdata", type=str,
                        dest="pathdata", help="Download data")
    parser.add_argument("-localrepo", "--localrepo", type=str,
                        dest="localrepo", help="Local repo")
    parser.add_argument("-pathout", "--pathout", type=str,
                        dest="pathout", help="Path output")
    args = parser.parse_args()

    if not args.pathdata:

        grb_input, redshift = download_data(year=args.grb_year)

    else:

        grb_input, redshift = use_downloaded_data(args.pathdata)

    for idx, val in enumerate(grb_input):

        dist_lum = cosm_calc(redshift[idx])

        grb_output = manipulate_txt(val)

        if grb_output == "":
            continue

        if args.condor:

            job_dir = args.pathdir + '/' + 'job_' + str(idx)

            os.mkdir(job_dir)

            write_bash_script_condor(
                grb_output, redshift[idx], dist_lum, job_dir, os.getenv('PWD'), args)

            write_file_sub(job_dir)

            subprocess.run(
                f"condor_submit -spool {job_dir}/sub", shell=True, check=True)

        else:

            write_bash_script(grb_output, redshift[idx], dist_lum, args)

            subprocess.run("bash script.sh", shell=True, check=True)


if __name__ == '__main__':
    main()
