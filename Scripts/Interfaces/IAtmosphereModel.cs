using UnityEngine;

namespace PhysicalSky
{
    public interface IAtmosphereModel
    {
        bool HasComputedData();
        bool NeedsRecompute();
        Color GetSunIlluminance(Vector3 positionRelativePlanetCenter, Vector3 sunDirection);
        void SetShaderUniforms(Material m);
        bool Compute(bool force = false);
        void ReleaseResources();
    }
}