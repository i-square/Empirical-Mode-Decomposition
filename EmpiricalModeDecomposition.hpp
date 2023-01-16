#pragma once

/*
 This code implements empirical mode decomposition in C.
 Required paramters include:
 - order: the number of IMFs to return
 - iterations: the number of iterations per IMF
 - locality: in samples, the nearest two extrema may be
 If it is not specified, there is no limit (locality = 0).

 Typical use consists of calling emdCreate(), followed by
 emdDecompose(), and then using the struct's "imfs" field
 to retrieve the data. Call emdClear() to deallocate memory
 inside the struct.
*/

#include <vector>

struct EMDData {
    int iterations = 1000;
    int order = 100;
    int locality = 0;
    float imf_thresh = 1e-8;
    int min_size = 0;
    int max_size = 0;
    std::vector<int> min_points;
    std::vector<int> max_points;
    std::vector<float> min;
    std::vector<float> max;
    std::vector<float> residual;

    void clear() {
        min_size = 0;
        max_size = 0;
        std::vector<int>().swap(min_points);
        std::vector<int>().swap(max_points);
        std::vector<float>().swap(min);
        std::vector<float>().swap(max);
        std::vector<float>().swap(residual);
    }
};

class EMD {
public:
    EMD();
    ~EMD();

public:
    int setup(
            int max_order = 100,
            int max_iterations = 1000,
            int locality = 0,
            float imf_thresh = 1e-8);
    std::vector<std::vector<float>> decompose(const std::vector<float>& signal);

private:
    inline int mirror_index(int i, int size);
    inline float abs_sum(const std::vector<float>& vec);
    bool check_peaks(const std::vector<float>& cur_imf);

    void make_extrema(const std::vector<float>& cur_imf);
    void interpolate(
            const std::vector<float>& in,
            std::vector<float>& out,
            std::vector<int>& points,
            int points_size);
    void update_imf(std::vector<float>& imf);
    void make_residue(const std::vector<float>& cur);
    void emd_clear();

private:
    EMDData* _emd;
};
