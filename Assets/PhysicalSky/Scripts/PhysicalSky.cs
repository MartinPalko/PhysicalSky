using UnityEngine;
using System;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    public class PhysicalSky : MonoBehaviour, IPhysicalSky
    {
        [SerializeField]
        Vector3 sunDirection = Vector3.one;
        public Vector3 SunDirection { get { return sunDirection; } set { sunDirection = value; } }

        [SerializeField]
        private float altitude = 0.1f;
        public float Altitude { get { return altitude; } set { altitude = value; } }

        [SerializeField]
        private AtmosphereModel atmosphere;
        public IAtmosphereModel Atmosphere { get { return atmosphere; } set { atmosphere = value as AtmosphereModel; } }

        [SerializeField]
        private Shader skyShader;

        private Material skyMaterial;

        private void Start()
        {
            skyMaterial = new Material(skyShader);
        }

        private void OnWillRenderObject()
        {
            Debug.Log("Will render object with camera " + Camera.current.name);
        }

        private void OnRenderObject()
        {
            if(Camera.current.name == "Reflection Probes Camera")
            {
                Debug.Log("Rendering reflection probes");
            }
        }

        private void Update()
        {
            if (atmosphere.NeedsRecompute)
                atmosphere.Compute();

            if (skyMaterial)
            {
                atmosphere.SetShaderUniforms(skyMaterial);
                skyMaterial.SetTexture("transmittance_texture", atmosphere.TransmittanceLUT);
                skyMaterial.SetTexture("scattering_texture", atmosphere.ScatteringLUT);
                skyMaterial.SetTexture("irradiance_texture", atmosphere.IrradianceLUT);

                skyMaterial.SetVector("camera", new Vector3(0, (atmosphere.PlanetaryRadius / 1000) + altitude, 0));

                RenderSettings.skybox = skyMaterial;
            }
        }
    }
}