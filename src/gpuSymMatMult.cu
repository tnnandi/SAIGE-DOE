#include <iostream>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"

#include "gpuSymMatMult.hpp"

#define IDX2C(i,j,ld)  (((j)*(ld))+(i))

gpuSymMatMult::gpuSymMatMult()
{
    m_isLoaded = false;
    m_rows = m_cols = 0;
    m_handle = nullptr;
    m_A = nullptr;
    m_x = nullptr;
    m_y = nullptr;
    m_z = nullptr;
}

gpuSymMatMult::~gpuSymMatMult()
{
    if (m_A) cudaFree(m_A);
    if (m_x) cudaFree(m_x);
    if (m_y) cudaFree(m_y);
    if (m_z) cudaFree(m_z);
    cublasDestroy((cublasHandle_t)m_handle);
}

int gpuSymMatMult::set_matrix(int rank, int n_rows, int n_cols, const float *A)
{
    cudaError_t cudaStat;
    cublasStatus_t cublasStat;
    size_t gpu_free, gpu_total;
    int device_id = -1;
    int ret = EXIT_SUCCESS;

    m_rows = n_rows;
    m_cols = n_cols;
    std::cout << "[" << rank << "] "
              << "A in nvcc = " << A << std::endl;

/*
    cudaStat = cudaSetDevice(rank);
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaSetDevice() has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
*/
    cudaStat = cudaGetDevice(&device_id);
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaGetDevice() has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    std::cout << "[" << rank << "] "
              << "GPU ID: " << device_id << std::endl;

    cudaStat = cudaMemGetInfo(&gpu_free, &gpu_total);
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMemGetInfo() has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    std::cout << "[" << rank << "] "
              << "GPU memory: " << gpu_free << " / " << gpu_total << std::endl;

    cudaStat = cudaMalloc((void**)&m_A, m_rows*m_cols*sizeof(*m_A));
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMalloc() for A has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    cudaStat = cudaMalloc((void**)&m_x, m_rows*sizeof(*m_x));
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMalloc() for x has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    cudaStat = cudaMalloc((void**)&m_y, m_cols*sizeof(*m_y));
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMalloc() for y has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    cudaStat = cudaMalloc((void**)&m_z, m_rows*sizeof(*m_z));
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMalloc() for z has failed: " << cudaStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    cublasStat = cublasCreate((cublasHandle_t *)&m_handle);
    if (cublasStat != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: cublasCreate() has failed: " << cublasStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }
    cublasStat = cublasSetMatrix(m_rows, m_cols, sizeof(*m_A), A, m_rows, m_A, m_rows);
    if (cublasStat != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: cublasSetMatrix() has failed: " << cublasStat << std::endl;
        ret = EXIT_FAILURE;
        goto out;
    }

    m_isLoaded = true;
out:
    if (ret == EXIT_FAILURE) {
        m_rows = m_cols = 0;
        if (m_A) cudaFree(m_A);
        if (m_x) cudaFree(m_x);
        if (m_y) cudaFree(m_y);
        if (m_z) cudaFree(m_z);
        if (m_handle) cublasDestroy((cublasHandle_t)m_handle);

        m_isLoaded = false;
        m_A = nullptr;
        m_x = nullptr;
        m_y = nullptr;
        m_z = nullptr;
        m_handle = nullptr;
    }
    return ret;
}

// Perform AA^T*x in two steps.
//   i)  y = A^T*x
//   ii) z = Ay
int gpuSymMatMult::sym_sgemv(int rank, int n_elem, const float *x, float *ret)
{
    cudaError_t cudaStat;
    cublasStatus_t cublasStat;
    float alpha = 1.0, beta = 0.0;

    if (n_elem != m_rows) {
        std::cerr << "[" << rank << "] "
                  << "Error: size mismatch: "
                  << "m_rows: " << m_rows << " n_elem: " << n_elem << std::endl;
        return EXIT_FAILURE;
    }

    memset(ret, 0, n_elem*sizeof(*ret));
    cudaStat = cudaMemcpy(m_x, x, n_elem*sizeof(*x), cudaMemcpyHostToDevice);
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMemcpy (host to Device) failed: " << cudaStat
                  << " m_x = " << m_x << " x = " << x << std::endl;
        return EXIT_FAILURE;
    }
    cublasStat = cublasSgemv((cublasHandle_t)m_handle, CUBLAS_OP_T, m_rows, m_cols, &alpha, m_A, m_rows, m_x, 1, &beta, m_y, 1);
    if (cublasStat != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: cublasSgemv failed (CUBLAS_OP_T): " << cublasStat << std::endl;
        return EXIT_FAILURE;
    }
    cublasStat = cublasSgemv((cublasHandle_t)m_handle, CUBLAS_OP_N, m_rows, m_cols, &alpha, m_A, m_rows, m_y, 1, &beta, m_z, 1);
    if (cublasStat != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "[" << rank << "] "
                  << "Error: cublasSgemv failed (CUBLAS_OP_N): " << cublasStat << std::endl;
        return EXIT_FAILURE;
    }
    cudaStat = cudaMemcpy(ret, m_z, n_elem*sizeof(*ret), cudaMemcpyDeviceToHost);
    if (cudaStat != cudaSuccess) {
        std::cerr << "[" << rank << "] "
                  << "Error: cudaMemcpy (device to host) failed: " << cudaStat << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
