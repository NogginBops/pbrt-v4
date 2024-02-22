
#ifndef PBRT_UTIL_REFLECTANCE_H
#define PBRT_UTIL_REFLECTANCE_H

#include <pbrt/pbrt.h>

#include <pbrt/util/taggedptr.h>

#include <pbrt/util/spectrum.h>

namespace pbrt {

class SpectrumReflectance;
class FluorescentReflectance;
class SampledReflectance;

class Reflectance : public TaggedPointer<SpectrumReflectance, FluorescentReflectance> {
  public:
	using TaggedPointer::TaggedPointer;
	std::string ToString() const;

	PBRT_CPU_GPU
    Float operator()(Float lambdaIn, Float lambdaOut) const;

    PBRT_CPU_GPU
    SampledReflectance Sample(const SampledWavelengths &lambda) const;
};

class SampledReflectance {
  public:
    PBRT_CPU_GPU
    SampledReflectance() {}
    // FIXME: Make this parameterized on NSpectrumSamples
    PBRT_CPU_GPU
    SampledReflectance(float f) { values == SquareMatrix<NSpectrumSamples>::Diag(f, f, f, f); }
    PBRT_CPU_GPU
    explicit SampledReflectance(SampledSpectrum f) {
        values = SquareMatrix<NSpectrumSamples>::Diag(f[0], f[1], f[2], f[3]);
    }
    PBRT_CPU_GPU
    SampledReflectance(SquareMatrix<NSpectrumSamples> m) { values = m; }

    PBRT_CPU_GPU
    static SampledReflectance FromSpectrum(const SampledSpectrum& s) {
        return SampledReflectance(s);
    }

    PBRT_CPU_GPU
    static SampledReflectance FromFloat(const Float& f) {
        return SampledReflectance(f);
    }

    PBRT_CPU_GPU
    static SampledSpectrum ToSpectrum(const SampledReflectance& r)
    {
        SampledSpectrum s;
#if 0
        // Approximate the albedo of a potentially fluorescent surface?
        for (int i = 0; i < NSpectrumSamples; i++) {
            for (int j = 0; j < NSpectrumSamples; j++) {
                s[i] += r.values[i][j];
            }
        }
#else
        // Take the diagonal of the reflectance matrix.
        for (int i = 0; i < NSpectrumSamples; i++) {
            s[i] = r.values[i][i];
        }
#endif
        return s;
    }

    PBRT_CPU_GPU
    SampledReflectance operator+(const SampledReflectance& r) const {
        return values + r.values;
    }

    PBRT_CPU_GPU
    SampledReflectance operator+=(const SampledReflectance& r) {
        values = values + r.values;
        return *this;
    }

    PBRT_CPU_GPU
    SampledReflectance operator+(const SampledSpectrum& r) {
        SampledReflectance ret = *this;
        for (int i = 0; i < NSpectrumSamples; i++) {
            ret.values[i][i] = ret.values[i][i] + r[i];
        }
        return ret;
    }

    PBRT_CPU_GPU
    SampledReflectance operator+=(const SampledSpectrum& r) {
        for (int i = 0; i < NSpectrumSamples; i++) {
            values[i][i] = values[i][i] + r[i];
        }
        return *this;
    }

    PBRT_CPU_GPU
    SampledReflectance operator*(const SampledReflectance& r) const {
        return values * r.values;
    }

    PBRT_CPU_GPU
    SampledReflectance operator*=(const SampledReflectance& r) {
        // FIXME: Multiplication order!
        values = values * r.values;
        return *this;
    }

    PBRT_CPU_GPU
    SampledSpectrum operator*(const SampledSpectrum& r) {
        SampledReflectance ret = *this;
        return ret *= r;
    }

    PBRT_CPU_GPU
    SampledSpectrum operator*=(const SampledSpectrum& r) {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; i++) {
            for (int j = 0; j < NSpectrumSamples; j++) {
                ret[i] += values[i][j] * r[j];
            }
        }
        return ret;
    }

    PBRT_CPU_GPU
    SampledReflectance MulDiag(const SampledSpectrum& r) {
        SampledReflectance ret = *this;
        for (int i = 0; i < NSpectrumSamples; i++) {
            for (int j = 0; j < NSpectrumSamples; j++) {
                ret.values[i][j] *= r[j];
            }
        }
        return ret;
    }

    PBRT_CPU_GPU
    SampledReflectance operator*(const Float f) const { return values * f; }

    PBRT_CPU_GPU
    SampledReflectance operator/(const Float f) { return values / f; }

    PBRT_CPU_GPU
    explicit operator bool() const {
        for (int i = 0; i < NSpectrumSamples; i++)
            for (int j = 0; j < NSpectrumSamples; j++)
                if (values[i][j] != 0)
                    return true;
        return false;
    }

    PBRT_CPU_GPU
    Float MaxComponentValue() const
    {
        Float max = -Infinity;
        for (int i = 0; i < NSpectrumSamples; i++) {
            for (int j = 0; j < NSpectrumSamples; j++) {
                max = std::max(max, values[i][j]);
            }
        }
        return max;
    }

    // FIXME: is operator() a good idea? Will it be confusing?
    PBRT_CPU_GPU
    Float operator()(int i, int j) const {
        DCHECK(i >= 0 && i < NSpectrumSamples);
        DCHECK(j >= 0 && j < NSpectrumSamples);
        return values[i][j];
    }

    PBRT_CPU_GPU
    Float& operator()(int i, int j) {
        DCHECK(i >= 0 && i < NSpectrumSamples);
        DCHECK(j >= 0 && j < NSpectrumSamples);
        return values[i][j];
    }

    std::string ToString() const;

  private:
    SquareMatrix<NSpectrumSamples> values;
};

PBRT_CPU_GPU inline SampledReflectance operator*(Float s, const SampledReflectance& m) {
    return m * s;
}

class SpectrumReflectance {
  public:
      // FIXME: Do we need to normalize 
    PBRT_CPU_GPU
    Float operator()(Float lambdaIn, Float lambdaOut) const { return lambdaIn == lambdaOut ? spectrum(lambdaIn) : 0.0f; }

    PBRT_CPU_GPU
    SampledReflectance Sample(const SampledWavelengths& lambda) const {
        SampledReflectance ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret(i, i) = spectrum(lambda[i]);
        return ret;
    }

    std::string ToString() const;

  private:
    Spectrum spectrum;
};

class FluorescentReflectance {
  public:
    PBRT_CPU_GPU
    Float operator()(Float lambdaIn, Float lambdaOut) const {
        // For now we will do bilinear interpolation of the reflectance matrix
        int inOffset = std::lround(lambdaIn) - lambda_min;
        int outOffset = std::lround(lambdaOut) - lambda_min;

        if (inOffset < 0 || inOffset >= reemission.size()) {
            return 0;
        }

        if (outOffset < 0 || outOffset >= reemission[inOffset].size())
        {
            return 0;
        }

        return reemission[inOffset][outOffset];
    }

    PBRT_CPU_GPU
    SampledReflectance Sample(const SampledWavelengths& lambda) const
    {
        SampledReflectance reflectance;
        // FIXME: Is this sampling even reasonable?
        for (int i = 0; i < NSpectrumSamples; i++) {
            int inOffset = std::lround(lambda[i]) - lambda_min;
            for (int j = 0; j < NSpectrumSamples; j++) {
                int outOffset = std::lround(lambda[i]) - lambda_min;
                if (inOffset < 0 || inOffset >= reemission.size() ||
                    outOffset < 0 || outOffset >= reemission[inOffset].size()) {
                    reflectance(i, j) = 0;
                } else {
                    reflectance(i, j) = reemission[inOffset][outOffset];
                }
            }
        }
        return reflectance;
    }

    std::string ToString() const;

  private:
    // FIXME: Fluorescence matrix data?
    int lambda_min, lambda_max;
    // FIXME: Do I need to mark this as GPU data for reemission[inOffset].size() to be
    // a __host__ __device__ function.
    std::vector<std::vector<Float>> reemission;
};

inline Float Reflectance::operator()(Float lambdaIn, Float lambdaOut) const {
    auto op = [&](auto ptr) { return (*ptr)(lambdaIn, lambdaOut); };
    return Dispatch(op);
}

inline SampledReflectance Reflectance::Sample(const SampledWavelengths &lambda) const {
    auto op = [&](auto ptr) { return ptr->Sample(lambda); };
    return Dispatch(op);
}

} // namespace pbrt

#endif  // PBRT_UTIL_REFLECTANCE_H