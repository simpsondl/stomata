#include <gpu_engine.hpp>
#include <stdexcept>

// CPU-only builds never dispatch here: the pipeline checks gpu_available().
// Direct callers receive an explicit error instead of misleading empty hits.
bool gpu_available() { return false; }

GpuMysersResult myers_gpu(const std::string&, GenomeView, size_t, size_t) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
GpuShiftAddResult shift_add_gpu(const std::string&, GenomeView, size_t, size_t) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
std::vector<GpuMysersResult> myers_gpu_batch_threshold(
    const std::vector<std::string>&, GenomeView, size_t, size_t, int) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
std::vector<GpuShiftAddResult> shift_add_gpu_batch_threshold(
    const std::vector<std::string>&, GenomeView, size_t, size_t, int) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
std::vector<GpuSparseResult> myers_gpu_batch_sparse(
    const std::vector<std::string>&, GenomeView, size_t, size_t, int) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
std::vector<GpuSparseResult> shift_add_gpu_batch_sparse(
    const std::vector<std::string>&, GenomeView, size_t, size_t, int) {
    throw std::runtime_error("Stomata was built without CUDA support");
}
