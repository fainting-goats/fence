
if [ $1 ]
then
    log=$1
else
    log="qexe.log"
fi

tail -n 30 $log
njobs=$(( $(qstat | wc -l)-2 ))
(($njobs<0)) && njobs=0
echo Remaining jobs: $njobs
qstat

# watch watchlog
#'tail -n 30 qexe.log; echo Remaining jobs: `qstat | wc -l`; qstat;'
