using UnityEngine;

namespace PhysicalSky
{
    public interface IAtmosphereModel
    {
        bool NeedsRecompute { get; }        
        void SetShaderUniforms(Material m);
        void Compute();
        void ReleaseResources();
    }
}