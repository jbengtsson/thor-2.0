# set current directory to working directory
#$ -cwd
# combine PBS standard output- and error files
#$ -j y

~/Thor-2.0/thor/wrk/nl_opt flat_file.dat 3
