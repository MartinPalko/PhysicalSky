/**
* Copyright (c) 2017 Eric Bruneton
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*
* Precomputed Atmospheric Scattering
* Copyright (c) 2008 INRIA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*/

using UnityEngine;

namespace PhysicalSky
{
    /// <summary>
    /// A 1:1 implementation of select features from AtmosphereLib.cginc on the CPU
    /// </summary>
    public static class AtmosphereLib
    {
        public struct DensityProfileLayer
        {
            public float width;
            public float exp_term;
            public float exp_scale;
            public float linear_term;
            public float constant_term;

            public DensityProfileLayer(float width, float exp_term, float exp_scale, float linear_term, float constant_term)
            {
                this.width = width;
                this.exp_term = exp_term;
                this.exp_scale = exp_scale;
                this.linear_term = linear_term;
                this.constant_term = constant_term;
            }

            public static bool operator ==(DensityProfileLayer x, DensityProfileLayer y)
            {
                return x.width == y.width &&
                    x.exp_term == y.exp_term &&
                    x.exp_scale == y.exp_scale &&
                    x.linear_term == y.linear_term &&
                    x.constant_term == y.constant_term;
            }
            public static bool operator !=(DensityProfileLayer x, DensityProfileLayer y) { return !x.Equals(y); }
            public override bool Equals(object obj) { return obj is DensityProfileLayer && this == (DensityProfileLayer)obj; }
            public override int GetHashCode() { return width.GetHashCode() ^ exp_term.GetHashCode() ^ exp_scale.GetHashCode() ^ linear_term.GetHashCode() ^ constant_term.GetHashCode(); }
        };

        public struct DensityProfile
        {
            public DensityProfileLayer layer0;
            public DensityProfileLayer layer1;

            public DensityProfile(DensityProfileLayer layer0, DensityProfileLayer layer1)
            {
                this.layer0 = layer0;
                this.layer1 = layer1;
            }

            public DensityProfile(DensityProfileLayer layer)
            {
                this.layer0 = layer;
                this.layer1 = layer;
            }

            public static bool operator ==(DensityProfile x, DensityProfile y)
            {
                return x.layer0 == y.layer0 && x.layer1 == y.layer1;
            }
            public static bool operator !=(DensityProfile x, DensityProfile y) { return !x.Equals(y); }
            public override bool Equals(object obj) { return obj is DensityProfile && this == (DensityProfile)obj; }
            public override int GetHashCode() { return layer0.GetHashCode() ^ layer1.GetHashCode(); }
        }

        [System.Serializable]
        public struct AtmosphereRenderParams
        {
            public Vector3 sky_spectral_radiance_to_luminance;
            public Vector3 sun_spectral_radiance_to_luminance;

            // The solar irradiance at the top of the atmosphere.
            public Vector3 solar_irradiance;
            // The sun's angular radius. Warning: the implementation uses approximations
            // that are valid only if this angle is smaller than 0.1 radians.
            public float sun_angular_radius;
            // The distance between the planet center and the bottom of the atmosphere.
            public float bottom_radius;
            // The distance between the planet center and the top of the atmosphere.
            public float top_radius;
            // The density profile of air molecules, i.e. a function from altitude to
            // dimensionless values between 0 (null density) and 1 (maximum density).
            public DensityProfile rayleigh_density;
            // The scattering coefficient of air molecules at the altitude where their
            // density is maximum (usually the bottom of the atmosphere), as a function of
            // wavelength. The scattering coefficient at altitude h is equal to
            // 'rayleigh_scattering' times 'rayleigh_density' at this altitude.
            public Vector3 rayleigh_scattering;
            // The density profile of aerosols, i.e. a function from altitude to
            // dimensionless values between 0 (null density) and 1 (maximum density).
            public DensityProfile mie_density;
            // The scattering coefficient of aerosols at the altitude where their density
            // is maximum (usually the bottom of the atmosphere), as a function of
            // wavelength. The scattering coefficient at altitude h is equal to
            // 'mie_scattering' times 'mie_density' at this altitude.
            public Vector3 mie_scattering;
            // The extinction coefficient of aerosols at the altitude where their density
            // is maximum (usually the bottom of the atmosphere), as a function of
            // wavelength. The extinction coefficient at altitude h is equal to
            // 'mie_extinction' times 'mie_density' at this altitude.
            public Vector3 mie_extinction;
            // The asymetry parameter for the Cornette-Shanks phase function for the
            // aerosols.
            public float mie_phase_function_g;
            // The density profile of air molecules that absorb light (e.g. ozone), i.e.
            // a function from altitude to dimensionless values between 0 (null density)
            // and 1 (maximum density).
            public DensityProfile absorption_density;
            // The extinction coefficient of molecules that absorb light (e.g. ozone) at
            // the altitude where their density is maximum, as a function of wavelength.
            // The extinction coefficient at altitude h is equal to
            // 'absorption_extinction' times 'absorption_density' at this altitude.
            public Vector3 absorption_extinction;
            // The average albedo of the ground.
            public Vector3 ground_albedo;
            // The cosine of the maximum Sun zenith angle for which atmospheric scattering
            // must be precomputed (for maximum precision, use the smallest Sun zenith
            // angle yielding negligible sky light radiance values. For instance, for the
            // Earth case, 102 degrees is a good choice - yielding mu_s_min = -0.2).
            public float mu_s_min;
        }

        private static float smoothstep(float edge0, float edge1, float x)
        {
            float t;
            t = Mathf.Clamp01((x - edge0) / (edge1 - edge0));
            return t * t * (3.0f - 2.0f * t);
        }

        private static float safeSqrt(float a)
        {
            return Mathf.Sqrt(Mathf.Max(a, 0.0f));
        }

        private static Vector3 SampleTex2D(Texture2D tex2D, Vector2 uv)
        {
            Color sample = tex2D.GetPixelBilinear(uv.x, uv.y);
            return new Vector3(sample.r, sample.g, sample.b);
        }

        private static Vector3 GetTransmittanceToTopAtmosphereBoundary(AtmosphereRenderParams atmosphere, Texture2D transmittanceLUT, float r, float mu)
        {
            Vector2 uv = GetTransmittanceTextureUvFromRMu(atmosphere, r, mu);
            return SampleTex2D(transmittanceLUT, uv);
        }

        private static float DistanceToTopAtmosphereBoundary(AtmosphereRenderParams atmosphere, float r, float mu)
        {
            Debug.Assert(r <= atmosphere.top_radius);
            Debug.Assert(mu >= -1.0f && mu <= 1.0f);
            float discriminant = r * r * (mu * mu - 1.0f) + atmosphere.top_radius * atmosphere.top_radius;
            return Mathf.Max(-r * mu + safeSqrt(discriminant), 0);
        }

        private static Vector2 GetTransmittanceTextureUvFromRMu(AtmosphereRenderParams atmosphere, float r, float mu)
        {
            Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
            Debug.Assert(mu >= -1.0 && mu <= 1.0);

            // Distance to top atmosphere boundary for a horizontal ray at ground level.
            float H = Mathf.Sqrt((atmosphere.top_radius * atmosphere.top_radius) - (atmosphere.bottom_radius * atmosphere.bottom_radius));
            // Distance to the horizon.
            float rho = safeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
            // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
            // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
            float d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
            float d_min = atmosphere.top_radius - r;
            float d_max = rho + H;
            float x_mu = (d - d_min) / (d_max - d_min);
            float x_r = rho / H;
            return new Vector2(x_mu, x_r);
        }

        private static Vector3 GetTransmittanceToSun(AtmosphereRenderParams atmosphere, Texture2D transmittanceLUT, float r, float mu_s)
        {
            float sin_theta_h = atmosphere.bottom_radius;
            float cos_theta_h = -Mathf.Sqrt(Mathf.Max(1.0f - sin_theta_h * sin_theta_h, 0.0f));
            return GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittanceLUT, r, mu_s) *
                smoothstep(-sin_theta_h * atmosphere.sun_angular_radius,
                    sin_theta_h * atmosphere.sun_angular_radius,
                    mu_s - cos_theta_h);
        }

        private static Vector3 GetSunIrradiance(AtmosphereRenderParams atmosphere, Texture2D transmittanceLUT, Vector3 positionRelativePlanetCenter, Vector3 sunDirection)
        {
            float r = positionRelativePlanetCenter.magnitude;
            float rmu = Vector3.Dot(positionRelativePlanetCenter, sunDirection);
            float mu_s = rmu / r;

            Vector3 transmittance = GetTransmittanceToSun(atmosphere, transmittanceLUT, r, mu_s);
            Vector3 irradiance = atmosphere.solar_irradiance;
            float f = smoothstep(-atmosphere.sun_angular_radius, atmosphere.sun_angular_radius, mu_s);

            return Vector3.Scale(irradiance, transmittance) * f;
        }

        public static Color GetSunIlluminance(AtmosphereRenderParams atmosphere, Texture2D transmittanceLUT, Vector3 positionRelativePlanetCenter, Vector3 sunDirection)
        {
            Vector3 sunIrradiance = GetSunIrradiance(atmosphere, transmittanceLUT, positionRelativePlanetCenter, sunDirection);
            Vector3 radianceToLuminance = atmosphere.sun_spectral_radiance_to_luminance;
            Vector3 sunIlluminance = Vector3.Scale(sunIrradiance, radianceToLuminance);

            return new Color(sunIlluminance.x, sunIlluminance.y, sunIlluminance.z, 1.0f);
        }
    }
}
