using UnityEngine;

namespace PhysicalSky
{
    public class SunControllerLightDriven : SunController
    {
        [SerializeField]
        protected Light sunDirectionalLight;

        protected void Update()
        {
            if (sky != null && sunDirectionalLight != null)
                sky.SunDirection = -sunDirectionalLight.transform.forward;
        }
    }
}