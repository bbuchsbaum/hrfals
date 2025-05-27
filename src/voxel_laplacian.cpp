#include <RcppArmadillo.h>
#include <array>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List build_voxel_laplacian_cpp(Rcpp::NumericVector volume, int connectivity = 6) {
    IntegerVector dims = volume.attr("dim");
    if(dims.size() != 3) {
        stop("'volume' must be 3D");
    }
    int nx = dims[0], ny = dims[1], nz = dims[2];
    int total = nx * ny * nz;
    const double* vol_ptr = volume.begin();

    std::vector<int> xs; xs.reserve(total);
    std::vector<int> ys; ys.reserve(total);
    std::vector<int> zs; zs.reserve(total);
    std::vector<int> index_map(total, 0);
    int idx = 0;
    for(int z = 0; z < nz; ++z) {
        for(int y = 0; y < ny; ++y) {
            for(int x = 0; x < nx; ++x) {
                int pos = x + nx * y + nx * ny * z;
                if(vol_ptr[pos] != 0) {
                    xs.push_back(x);
                    ys.push_back(y);
                    zs.push_back(z);
                    index_map[pos] = ++idx; // 1-based
                }
            }
        }
    }

    int v = idx;
    if(v == 0) {
        stop("mask contains no voxels");
    }

    std::vector<std::array<int,3>> offsets;
    for(int dx = -1; dx <= 1; ++dx) {
        for(int dy = -1; dy <= 1; ++dy) {
            for(int dz = -1; dz <= 1; ++dz) {
                if(dx == 0 && dy == 0 && dz == 0) continue;
                int sumabs = std::abs(dx) + std::abs(dy) + std::abs(dz);
                if((connectivity == 6 && sumabs == 1) ||
                   (connectivity == 18 && (sumabs == 1 || sumabs == 2)) ||
                   (connectivity == 26)) {
                    offsets.push_back({{dx, dy, dz}});
                }
            }
        }
    }

    std::vector<arma::uword> I;
    std::vector<arma::uword> J;
    I.reserve(offsets.size() * v);
    J.reserve(offsets.size() * v);

    for(const auto& off : offsets) {
        int dx = off[0];
        int dy = off[1];
        int dz = off[2];
        for(int k = 0; k < v; ++k) {
            int x = xs[k] + dx;
            int y = ys[k] + dy;
            int z = zs[k] + dz;
            if(x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) continue;
            int pos = x + nx * y + nx * ny * z;
            int nb = index_map[pos];
            if(nb > 0 && (k + 1) < nb) {
                arma::uword from = k;
                arma::uword to = nb - 1;
                I.push_back(from);
                J.push_back(to);
                I.push_back(to);
                J.push_back(from);
            }
        }
    }

    arma::sp_mat Adj(v, v);
    if(!I.empty()) {
        std::vector<double> val(I.size(), 1.0);
        Adj = arma::sp_mat(I.begin(), J.begin(), val.begin(), v, v);
    }

    arma::vec degree = arma::sum(Adj, 1);
    arma::sp_mat L = -Adj;
    for(arma::uword i = 0; i < degree.n_elem; ++i) {
        L(i, i) = degree(i);
    }

    return Rcpp::List::create(_["L"] = L, _["degree"] = degree);
}

