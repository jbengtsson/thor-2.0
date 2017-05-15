# set current directory to working directory
#$ -cwd
# combine PBS standard output- and error files
#$ -j y

# sm3 with b_4.
~/Thor-2.0/thor/wrk/nl_cg flat_file.dat 1
