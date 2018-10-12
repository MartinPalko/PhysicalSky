using UnityEngine;

namespace PhysicalSky
{
    public interface IAtmosphereModel
    {
        bool NeedsRecompute { get; }
        void SetShaderUniforms(Material m);
        bool Compute(bool force = false);
        void ReleaseResources();
    }
}