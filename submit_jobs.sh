#!/bin/bash

for xtal in COF-103.cif SBMOF-1.cssr NiPyC2.cif
do
    echo "submitting job for $xtal"
    sbatch -J $xtal -A simon-grp -c4\
    -o "$xtal.o" -e "$xtal.e" --export=xtal="$xtal" gcmc_submit.sh
done

# NOTE: sbatch does not launch tasks, it requests an allocation of resources and submits a batch script
# optional arguements can be found at https://slurm.schedmd.com/sbatch.html

## SBATCH -J job_name # name that will apprear with id upon query
## SBATCH -A simoncor # sponsored acount = research group
## SBATCH -p mime5 # name of queue or partition to run on
## SBATCH -n 4 # (--ntasks) how many job steps run within the allocated resources
## SBATCH --cpus-per-task=1 # how many processors per task | default=1)
##     arguement is set in julia submission as the "-p" flag, so this can be left as default.
## SBATCH -o filename.out # output file
## SBATCH -e filename.err # error file
## SBATCH --mail-type=BEGIN,END,FAIL # event for mailing user
## SBATCH --mail-user=myONID@oregonstate.edu # who will receive the email for the specified events
