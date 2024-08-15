#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_GRID=`pwd`/../main
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
INPUTDIR1=/data/lihl/code/2d-grid-generate/hyperbolic/project1/output
INPUTDIR2=/data/lihl/code/2d-grid-generate/hyperbolic/project2/output
STRETCH_FILE1=`pwd`/arc_len_file1.txt
STRETCH_FILE2=`pwd`/arc_len_file2.txt

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p ${PROJDIR}
mkdir -p ${OUTPUT_DIR}

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "input_grids_info" : [
    {
      "grid_import_dir" : "${INPUTDIR1}",
      "number_of_grid_points" : [200,400],
      "number_of_mpiprocs_in" : [2,2],
      "flag_stretch" : 0,
      "stretch_file" : "${STRETCH_FILE1}"
    },
    {
      "grid_import_dir" : "${INPUTDIR2}",
      "number_of_grid_points" : [200,400],
      "number_of_mpiprocs_in" : [2,2],
      "flag_stretch" : 0,
      "stretch_file" : "${STRETCH_FILE2}"
    }
  ],
    
  "stretch_direction" : "x",
  "merge_direction" : "x",

  "number_of_mpiprocs_out" : [1,2],

  "check_orth" : 1,
  "check_jac" : 1,
  "check_ratio" : 1,
  "check_step_xi" : 1,
  "check_step_zt" : 1,
  "check_smooth_xi" : 1,
  "check_smooth_zt" : 1,

  "flag_sample" : 0,
  "sample_factor_xi" : 1,
  "sample_factor_zt" : 1,

  "grid_export_dir" : "${OUTPUT_DIR}"

  "flag_pml" : 0,
  "pml_layers" : {
         "number_of_pml_x1" : 10,
         "number_of_pml_x2" : 10,
         "number_of_pml_z1" : 10,
         "number_of_pml_z2" : 0
  }

}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- gen run script
cat << ieof > ${PROJDIR}/grid_generate.sh
#!/bin/bash

set -e

printf "\nStart grid post process ...\n";
time $EXEC_GRID $PAR_FILE 2>&1 |tee log
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
