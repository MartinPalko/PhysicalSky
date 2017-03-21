using UnityEngine;

namespace PhysicalSky
{
    public class SunControllerLightDriven : SunController
    {
        [SerializeField]
        protected Light sunDirectionalLight;

        protected void Update()
        {
            sky.SunDirection = -sunDirectionalLight.transform.forward;
        }
    }
}