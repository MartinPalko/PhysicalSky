using UnityEngine;

namespace PhysicalSky
{
    public interface IPhysicalSky
    {
        Vector3 SunDirection { get; set; }
        float Altitude { get; set; }
        IAtmosphereModel Atmosphere { get; set; }
    }
}