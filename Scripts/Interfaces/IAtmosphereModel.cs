using UnityEngine;

namespace PhysicalSky
{
    public interface IAtmosphereModel
    {
        bool HasComputedData();
        bool NeedsRecompute();
        void SetShaderUniforms(Material m);
        bool Compute(bool force = false);
        void ReleaseResources();
    }
}