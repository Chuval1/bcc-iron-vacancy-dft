#!/bin/bash
# Master script to submit all convergence test jobs

# Submit 1x1x1 with k-mesh k16x16x16
cd 1x1x1/k16x16x16
sbatch submit.sh fe_bulk_1x1x1_k161616.in fe_bulk_1x1x1_k161616.out
cd ../../

# Submit 1x1x1 with k-mesh k20x20x20
cd 1x1x1/k20x20x20
sbatch submit.sh fe_bulk_1x1x1_k202020.in fe_bulk_1x1x1_k202020.out
cd ../../

# Submit 1x1x1 with k-mesh k24x24x24
cd 1x1x1/k24x24x24
sbatch submit.sh fe_bulk_1x1x1_k242424.in fe_bulk_1x1x1_k242424.out
cd ../../

# Submit 2x2x2 with k-mesh k8x8x8
cd 2x2x2/k8x8x8
sbatch submit.sh fe_bulk_2x2x2_k888.in fe_bulk_2x2x2_k888.out
cd ../../

# Submit 2x2x2 with k-mesh k10x10x10
cd 2x2x2/k10x10x10
sbatch submit.sh fe_bulk_2x2x2_k101010.in fe_bulk_2x2x2_k101010.out
cd ../../

# Submit 2x2x2 with k-mesh k12x12x12
cd 2x2x2/k12x12x12
sbatch submit.sh fe_bulk_2x2x2_k121212.in fe_bulk_2x2x2_k121212.out
cd ../../

# Submit 3x3x3 with k-mesh k6x6x6
cd 3x3x3/k6x6x6
sbatch submit.sh fe_bulk_3x3x3_k666.in fe_bulk_3x3x3_k666.out
cd ../../

# Submit 3x3x3 with k-mesh k8x8x8
cd 3x3x3/k8x8x8
sbatch submit.sh fe_bulk_3x3x3_k888.in fe_bulk_3x3x3_k888.out
cd ../../

# Submit 3x3x3 with k-mesh k10x10x10
cd 3x3x3/k10x10x10
sbatch submit.sh fe_bulk_3x3x3_k101010.in fe_bulk_3x3x3_k101010.out
cd ../../

# Submit 4x4x4 with k-mesh k4x4x4
cd 4x4x4/k4x4x4
sbatch submit.sh fe_bulk_4x4x4_k444.in fe_bulk_4x4x4_k444.out
cd ../../

# Submit 4x4x4 with k-mesh k6x6x6
cd 4x4x4/k6x6x6
sbatch submit.sh fe_bulk_4x4x4_k666.in fe_bulk_4x4x4_k666.out
cd ../../

# Submit 4x4x4 with k-mesh k8x8x8
cd 4x4x4/k8x8x8
sbatch submit.sh fe_bulk_4x4x4_k888.in fe_bulk_4x4x4_k888.out
cd ../../

echo "All jobs submitted."
