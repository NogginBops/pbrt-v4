

#include <pbrt/util/reflectance.h>

namespace pbrt {

std::string Reflectance::ToString() const {
    if (!ptr())
        return "(nullptr)";

    auto tostr = [&](auto ptr) { return ptr->ToString(); };
    return DispatchCPU(tostr);
}

std::string SpectrumReflectance::ToString() const {
    return StringPrintf("[ SpectrumReflectance spectrum: %s ]", spectrum.ToString());
}

std::string FluorescentReflectance::ToString() const {
    // FIXME: Print the matrix?
    return StringPrintf("[ FluorescentReflectance TODO ]");
}

std::string SampledReflectance::ToString() const {
    return StringPrintf("[ SampledReflectance matrix: %s ]", values.ToString());
}

}
