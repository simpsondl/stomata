#!/bin/bash
# GPU Safety Test Harness
# Run this before every GPU code execution

set -e

echo "=== GPU Safety Checks ==="

# 1. Check GPU availability
if ! nvidia-smi &> /dev/null; then
    echo "ERROR: nvidia-smi not available. GPU may not be accessible."
    exit 1
fi

# 2. Check GPU memory
FREE_MEM=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -n 1)
if [ "$FREE_MEM" -lt 4000 ]; then
    echo "WARNING: Less than 4GB free GPU memory ($FREE_MEM MB)"
    echo "Consider freeing GPU resources before continuing"
fi

# 3. Check for zombie GPU processes
GPU_PROCS=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
if [ "$GPU_PROCS" -gt 5 ]; then
    echo "WARNING: $GPU_PROCS GPU processes detected. May cause resource contention."
fi

# 4. Test CUDA compilation
echo "Testing CUDA compiler..."
cat > /tmp/cuda_test.cu << 'CUDA_EOF'
#include <cuda_runtime.h>
#include <stdio.h>

__global__ void test_kernel() { }

int main() { 
    test_kernel<<<1,1>>>(); 
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
        return 1;
    }
    printf("CUDA test successful\n");
    return 0; 
}
CUDA_EOF

nvcc -arch=sm_89 /tmp/cuda_test.cu -o /tmp/cuda_test 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: CUDA compilation failed"
    exit 1
fi

/tmp/cuda_test
if [ $? -ne 0 ]; then
    echo "ERROR: CUDA test kernel failed"
    exit 1
fi

rm -f /tmp/cuda_test /tmp/cuda_test.cu

echo "✓ All GPU safety checks passed"
echo "Free GPU Memory: ${FREE_MEM} MB"
