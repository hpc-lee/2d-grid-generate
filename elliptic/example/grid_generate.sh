#!/bin/bash

#set -x
set -e

date

#-- system related dir
MPIDIR=/data/apps/openmpi/4.1.5-cuda-aware
#-- program related dir
EXEC_GRID=`pwd`/../main
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
INPUTDIR=`pwd`

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p ${PROJDIR}
mkdir -p ${OUTPUT_DIR}

#-- total x mpi procs
NPROCS_X=1
#-- total z mpi procs
NPROCS_Z=1
#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "number_of_grid_points_x" : 1708,
  "number_of_grid_points_z" : 601,

  "number_of_mpiprocs_x" : $NPROCS_X,
  "number_of_mpiprocs_z" : $NPROCS_Z,

  "check_orth" : 1,
  "check_jac" : 1,
  "check_ratio" : 1,
  "check_step_xi" : 1,
  "check_step_zt" : 1,
  "check_smooth_xi" : 1,
  "check_smooth_zt" : 1,

  "geometry_input_file" : "${INPUTDIR}/data_file_2d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "grid_method" : {
      "#tfi":"",
      "#elli_diri" : {
          "coef" : [20,20,20,20],
          "iter_err" : 1E-2,
          "max_iter" : 5E3
      },
      "elli_higen" : {
          "coef" : [1000,1000,20,20],
          "iter_err" : 1E-2,
          "max_iter" : 5E3
      }
  }
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- get np
NUMPROCS_X=`grep number_of_mpiprocs_x ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Z=`grep number_of_mpiprocs_z ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS=$(( NUMPROCS_X*NUMPROCS_Z ))
echo $NUMPROCS_X $NUMPROCS_Z $NUMPROCS

#-- gen run script
cat << ieof > ${PROJDIR}/grid_generate.sh
#!/bin/bash

set -e

printf "\nUse $NUMPROCS CPUs on following nodes:\n"

printf "\nStart grid generate ...\n";
time $MPIDIR/bin/mpiexec -np $NUMPROCS $EXEC_GRID $PAR_FILE 100 2 2>&1 |tee log
if [ $? -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/grid_generate.sh
${PROJDIR}/grid_generate.sh

date

# vim:ts=4:sw=4:nu:et:ai:
