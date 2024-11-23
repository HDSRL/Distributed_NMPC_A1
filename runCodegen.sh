#!/bin/bash

# Navigate to the SRB_NMPC_Casadi_Formulation/build directory, create it if it does not exist
if [ ! -d "SRB_NMPC_Casadi_Formulation/build" ]; then
  mkdir -p SRB_NMPC_Casadi_Formulation/build
fi
cd SRB_NMPC_Casadi_Formulation/build

# Run make with 16 parallel jobs
make -j16

# Execute the nmpccodegen executable
./nmpccodegen

# Copy the specified files to the target directory
cp nlp_code_f.so ../../build/
cp nlp_code.so ../../build/

# Print a message indicating the script has completed
echo "Build and file copy completed successfully."
