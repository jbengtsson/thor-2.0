# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$HOME/git_repos/thor-2.0/projects/src/tbs $lat_file
