#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_GRID=`pwd`/../main_grid_2d
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

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{

  "check_orth" : 1,
  "check_jac" : 1,
  "check_ratio" : 1,
  "check_step" : 1,
  "check_smooth" : 1,

  "flag_strech_x" : 0,
  "strech_x_coef" : 0.0001,
  "flag_strech_z" : 1,
  "strech_z_coef" : 0.0001,

  "flag_sample_x" : 1,
  "sample_factor_x" : 1.5,
  "flag_sample_x" : 1,
  "sample_factor_z" : 1.5,

  "geometry_input_file" : "${INPUTDIR}/data_file_2d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "grid_method" : {
      "#linear_TFI" : "",
      "hermite" : {
          "coef" : 0.3
      },
      "#elli_diri" : {
          "coef" : 0.3
      },
      "#elli_higen" : {
          "coef" : 0.3
      },
      "#parabolic" : {
          "coef" : 0.3
      },
      "#hyperbolic" : {
          "coef" : 0.3,
          "num_layers" : 30,
          "step_input_file" : "${INPUTDIR}/step_file_2d.txt"
      }
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

printf "\nStart grid generate ...\n";
time $EXEC_GRID $PAR_FILE 100 2>&1 |tee log
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
