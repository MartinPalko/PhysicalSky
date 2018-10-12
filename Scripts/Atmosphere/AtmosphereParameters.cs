using UnityEngine;

namespace PhysicalSky
{
    [System.Serializable]
    public struct AtmosphereParameters
    {
        public enum LuminanceType
        {
            // Render the spectral radiance at kLambdaR, kLambdaG, kLambdaB.
            none = 0,
            // Render the sRGB luminance, using an approximate (on the fly) conversion
            // from 3 spectral radiance values only (see section 14.3 in "A Qualitative
            // and Quantitative Evaluation of 8 Clear Sky Models"
            // https://arxiv.org/pdf/1612.04336.pdf).
            approximate = 1,
            // Render the sRGB luminance, precomputed from 15 spectral radiance values
            // (see section 4.4 in "Real-time Spectral Scattering in Large-scale 
            // Natural Participating Media"
            // http://www.oskee.wz.cz/stranka/uploads/SCCG10ElekKmoch.pdf).
            precomputed = 2
        }

        public int scatteringOrders;
        public LuminanceType luminance;
        public bool useOzone;
        public bool useConstantSolarSpectrum;
        public float constantSolarIrradiance;
        public float sunAngularRadius;
        public float planetaryRadius;
        public float atmosphereThickness;
        public float rayleigh;
        public float rayleighScaleHeight;
        public float mieScaleHeight;
        public float mieAngstromAlpha;
        public float mieAngstromBeta;
        public float mieSingleScatteringAlbedo;
        public float miePhaseFunctionG;
        public float groundAlbedo;
        public float maxSunZenithAngle;

        public override int GetHashCode()
        {
            return
                scatteringOrders.GetHashCode() ^
                luminance.GetHashCode() ^
                useOzone.GetHashCode() ^
                useConstantSolarSpectrum.GetHashCode() ^
                constantSolarIrradiance.GetHashCode() ^
                sunAngularRadius.GetHashCode() ^
                planetaryRadius.GetHashCode() ^
                atmosphereThickness.GetHashCode() ^
                rayleigh.GetHashCode() ^
                rayleighScaleHeight.GetHashCode() ^
                mieScaleHeight.GetHashCode() ^
                mieAngstromAlpha.GetHashCode() ^
                mieAngstromBeta.GetHashCode() ^
                mieSingleScatteringAlbedo.GetHashCode() ^
                miePhaseFunctionG.GetHashCode() ^
                groundAlbedo.GetHashCode() ^
                maxSunZenithAngle.GetHashCode();
        }

        public static bool operator ==(AtmosphereParameters x, AtmosphereParameters y)
        {
            return 
                x.scatteringOrders == y.scatteringOrders &&
                x.luminance == y.luminance &&
                x.useOzone == y.useOzone &&
                x.useConstantSolarSpectrum == y.useConstantSolarSpectrum &&
                x.constantSolarIrradiance == y.constantSolarIrradiance &&
                x.sunAngularRadius == y.sunAngularRadius &&
                x.planetaryRadius == y.planetaryRadius &&
                x.atmosphereThickness == y.atmosphereThickness &&
                x.rayleigh == y.rayleigh &&
                x.rayleighScaleHeight == y.rayleighScaleHeight &&
                x.mieScaleHeight == y.mieScaleHeight &&
                x.mieAngstromAlpha == y.mieAngstromAlpha &&
                x.mieAngstromBeta == y.mieAngstromBeta &&
                x.mieSingleScatteringAlbedo == y.mieSingleScatteringAlbedo &&
                x.miePhaseFunctionG == y.miePhaseFunctionG &&
                x.groundAlbedo == y.groundAlbedo &&
                x.maxSunZenithAngle == y.maxSunZenithAngle;
        }

        public static bool operator !=(AtmosphereParameters x, AtmosphereParameters y) { return !x.Equals(y); }
        public override bool Equals(object obj) { return obj is AtmosphereParameters && this == (AtmosphereParameters)obj; }

        public static AtmosphereParameters Validate(AtmosphereParameters parameters)
        {
            parameters.scatteringOrders = Mathf.Clamp(parameters.scatteringOrders, 0, 10);
            parameters.constantSolarIrradiance = Mathf.Max(0.001f, parameters.constantSolarIrradiance);
            parameters.sunAngularRadius = Mathf.Clamp(parameters.sunAngularRadius, 0.001f, 0.1f);
            parameters.planetaryRadius = Mathf.Max(1.0f, parameters.planetaryRadius);
            parameters.atmosphereThickness = Mathf.Max(1.0f, parameters.atmosphereThickness);
            parameters.rayleigh = Mathf.Max(1.0e-22f, parameters.rayleigh);
            parameters.rayleighScaleHeight = Mathf.Max(0.0f, parameters.rayleighScaleHeight);
            parameters.mieScaleHeight = Mathf.Max(0.0f, parameters.mieScaleHeight);
            parameters.mieAngstromAlpha = Mathf.Max(0.0f, parameters.mieAngstromAlpha);
            parameters.mieAngstromBeta = Mathf.Max(0.0f, parameters.mieAngstromBeta);
            parameters.mieSingleScatteringAlbedo = Mathf.Max(0.0f, parameters.mieSingleScatteringAlbedo);
            parameters.miePhaseFunctionG = Mathf.Max(0.0f, parameters.miePhaseFunctionG);
            parameters.groundAlbedo = Mathf.Max(0.0f, parameters.groundAlbedo);
            parameters.maxSunZenithAngle = Mathf.Max(0.0f, parameters.maxSunZenithAngle);
            return parameters;
        }

        public static AtmosphereParameters defaultEarth
        {
            get
            {
                AtmosphereParameters earth;
                earth.scatteringOrders = 4;
                earth.luminance = LuminanceType.precomputed;
                earth.useOzone = true;
                earth.useConstantSolarSpectrum = false;
                earth.constantSolarIrradiance = 1.5f;
                earth.sunAngularRadius = 0.00872665f;
                earth.planetaryRadius = 6360000.0f;
                earth.atmosphereThickness = 60000.0f;
                earth.rayleigh = 1.24062e-6f;
                earth.rayleighScaleHeight = 8000.0f;
                earth.mieScaleHeight = 1200.0f;
                earth.mieAngstromAlpha = 0.0f;
                earth.mieAngstromBeta = 5.328e-3f;
                earth.mieSingleScatteringAlbedo = 0.9f;
                earth.miePhaseFunctionG = 0.8f;
                earth.groundAlbedo = 0.3f;
                earth.maxSunZenithAngle = 102.0f / 180.0f * Mathf.PI;
                return earth;
            }
        }
    }
}