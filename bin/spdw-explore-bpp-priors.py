#! /usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bpp_control_filepath",
            action="store",
            help="Path to BPP control file.")
    args = parser.parse_args()
    control_filepath = os.path.expanduser(os.path.expandvars(args.bpp_control_filepath))
    control_file_dir = os.path.dirname(os.path.abspath(control_filepath))
    with open(control_filepath) as src:
        control_file_contents = src.read()
    theta_prior_a = 3.0
    tau_prior_a = 3.0
    job_files = []
    job_root_dir = os.path.join(os.path.abspath(os.path.curdir), "analyses")
    os.makedirs(job_root_dir, exist_ok=True)
    seqfile_pattern =     re.compile(r"seqfile\s+=\s+(\S+)", re.MULTILINE | re.IGNORECASE)
    imapfile_pattern =    re.compile(r"Imapfile\s+=\s+(\S+)", re.MULTILINE | re.IGNORECASE)
    theta_prior_pattern = re.compile(r"thetaprior\s+=\s+([0-9.]+)\s+([0-9.]+)", re.MULTILINE | re.IGNORECASE)
    tau_prior_pattern =   re.compile(r"tauprior\s+=\s+([0-9.]+)\s+([0-9.]+)", re.MULTILINE | re.IGNORECASE)
    for theta_b_exp in range(1, 7):
        for tau_b_exp in range(1, 7):
            theta_prior_b = 2.0 / (10**theta_b_exp)
            theta_prior_mean = theta_prior_b / (theta_prior_a - 1)
            tau_prior_b = 2.0 / (10**tau_b_exp)
            tau_prior_mean = tau_prior_b / (tau_prior_a - 1)
            title = "theta_{:.0e}_tau_{:.0e}".format(theta_prior_b, tau_prior_b)
            # print("{}: thetaprior = {:e} {:e} (mean = {:e}); tauprior = {} {} (mean = {:e})".format(
            #         title,
            #         theta_prior_a,
            #         theta_prior_b,
            #         theta_prior_mean,
            #         tau_prior_a,
            #         tau_prior_b,
            #         tau_prior_mean))
            print(title)
            job_dir = os.path.join(job_root_dir, title)
            os.makedirs(job_dir, exist_ok=True)
            job_control = control_file_contents
            job_control = seqfile_pattern.sub(r"seqfile = {}/\g<1>".format(os.path.join(control_file_dir)), job_control)
            job_control = imapfile_pattern.sub(r"Imapfile = {}/\g<1>".format(os.path.join(control_file_dir)), job_control)
            job_control = theta_prior_pattern.sub(r"thetaprior = \g<1> {} * (mean = {})".format(theta_prior_b, theta_prior_mean), job_control)
            job_control = tau_prior_pattern.sub(r"tauprior = \g<1> {} * (mean = {})".format(tau_prior_b, tau_prior_mean), job_control)
            job_filename = os.path.join(job_dir, "bpp.ctl")
            with open(job_filename, "w") as dest:
                dest.write(job_control)
            sys.exit()

if __name__ == "__main__":
    main()
