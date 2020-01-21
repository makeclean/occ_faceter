#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#define USE_HALF_FLOAT 1
#if USE_HALF_FLOAT
#include "UniqueId/half.hpp"
#else
#error "in-house double to unit16_t mapping is not yet implemented"
#endif

/* *
 * these functions are used "GeometryPropertyBuilder" to calc shape unique ID
 * from 4 key values: volume (for solid) or surface area, plus centre of mass cooordiante (x, y, z)
 * all geometry (shapes) are saved as a compound of solids, for each solid readback,
 * using methods in "GeometryPropertyBuilder" to calc 4 key values. then calc unique ID
 * using this unique ID, to match and locate meta information (json file) for each solid
 *
 * dependencies: half.hpp, json.hpp, OCCT
 * */

namespace PPP
{
    namespace Utilities
    {
        // todo: typedef uint64_t UniqueIdType;

        /// small endianness is more popular in CPU worlds
        bool isBigEndian(void)
        {
            union {
                uint32_t i;
                char c[4];
            } bint = {0x01020304};

            return bint.c[0] == 1;
        }

        /// for unsigned int only
        template <class T> static T round(T value, T tol)
        {
            T remained = value % (tol);
            T r = value / (tol);
            if (remained >= (tol / 2))
                r += 1;
            return r * (tol);
        }

        /// map from double float to half precision float (binary16)
        /// with EPS approximately 0.001, so this mapping uses log10 not equal space mapping
        /// bit_mask can be used to further reduced EPS
        /// then from half float to underneath uint16_t by reinterpret_cast
        /// this function use third-party lib: HalfFloat,
        /// if the value is close to zero (zeroThreshold = 1e-3), then set it as zero,
        std::uint16_t double2uint16(double value)
        {
            static double zeroThreshold = 1e-4;
            std::uint16_t round_precision = 0x0008; // mask out extra significant bits
            if (std::fabs(value) < zeroThreshold)
                value = 0.0;
#if USE_HALF_FLOAT
            using namespace half_float;
            half hf(static_cast<float>(value));
            std::uint16_t* p = reinterpret_cast<std::uint16_t*>(&hf);
            std::uint16_t ret = round(*p, round_precision);
            return ret;
#endif
        }


        /// due to unlikely floating point error, calculated ID can be different after rounding
        std::vector<uint64_t> nearbyIds(uint64_t id)
        {
            std::vector<uint64_t> nids;
            // todo:
            return nids;
        }

        // approximate equal by float point data comparison
        bool uniqueIdEqual(uint64_t id1, std::uint64_t id2)
        {
            std::uint16_t round_precision = 0x0008; // defined in `double2uint16()`
            for (size_t i = 0; i < 4; i++)
            {
                uint16_t i1 = (id1 << i * 16) & 0xFFFF;
                uint16_t i2 = (id2 << i * 16) & 0xFFFF;
                // todo: this does not consider overflow
                if (i1 > i2 + round_precision || i1 < i2 - round_precision)
                {
                    return false;
                }
            }
            return true;
        }

        /// todo: endianess problem.
        /// this should generate unique ID for a vector of 4 double values
        /// this is not universal unique ID (UUID) yet
        std::uint64_t uniqueId(const std::vector<double> values)
        {
            std::uint64_t ret = 0x00000000;
            assert(values.size() == 4);

            for (size_t i = 0; i < values.size(); i++)
            {
                std::uint64_t tmp = double2uint16(values[i]);
                ret ^= (tmp << i * 16);
            }
            return ret;
        }

        /// from float array to int[] then into a hash value
        /// using cm3, cm as unit, assuming geometry has mm as unit
        std::uint64_t geometryUniqueID(double volume, std::vector<double> centerOfMass)
        {
            const double scale = 1e2; // should equal to that in GeometryPropertyBuilder class in parallel-preprocessor
            double gp = volume / (scale * scale * scale);
            std::vector<double> pv{gp, centerOfMass[0] / scale, centerOfMass[1] / scale, centerOfMass[2] / scale};
            return uniqueId(pv);
        }
    } // namespace Utilities
} // namespace PPP