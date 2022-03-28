#ifndef GPU_SYM_MAT_MULT_HPP
#define GPU_SYM_MAT_MULT_HPP

class gpuSymMatMult
{
public:
    bool m_isLoaded;
    size_t m_rows;
    size_t m_cols;
    size_t m_global_cols;
    void *m_handle;
    float *m_A;
    float *m_x;
    float *m_y;
    float *m_z;

    gpuSymMatMult();
    ~gpuSymMatMult();

    int set_matrix(int rank, size_t n_rows, size_t n_cols, const float *A);
    int sym_sgemv(int rank, size_t n_elem, const float *x, float *ret);
};
#endif // GPU_SYM_MAT_MULT_HPP
