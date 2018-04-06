# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

echo $ns

$HOME/git_repos/thor-2.0/projects/src/tune_conf_2 flat_file.dat $ns
