from download import download_data
from manipulate_txt import manipulate_txt
from cosm_calc import cosm_calc
import subprocess
from argparse import ArgumentParser


def write_bash_script(file, redshift, dist_lum):

    with open('script.sh', 'w') as script:
        script.write("#!/usr/bin/env bash\n")
        script.write("source activate JetFit\n")
        script.write(
            f"python fitter.py --grb {file} --z {redshift} --dL {dist_lum}")
    subprocess.run("chmod +x script.sh", shell=True, check=True)


def main(args=None):

    parser = ArgumentParser(description='Pipeline for GRB analysis')

    parser.add_argument("-grb", "--grb", type=str,
                        dest="grb_name", help="GRB name")
    parser.add_argument("-year", "--year", type=int,
                        dest="grb_year", help="GRB year")
    args = parser.parse_args()

    grb_input, redshift = download_data(year=args.grb_year)

    for idx, val in enumerate(grb_input):

        dist_lum = cosm_calc(redshift[idx])

        grb_output = manipulate_txt(val)

        if grb_output == "":
            continue

        write_bash_script(grb_output, redshift[idx], dist_lum)

        subprocess.run("bash script.sh", shell=True, check=True)


if __name__ == '__main__':
    main()
