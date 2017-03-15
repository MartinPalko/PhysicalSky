using System;
using UnityEngine;

public class PhysicalSky : MonoBehaviour
{
    public AtmosphereModel atmosphereModel;
    
    private static Mesh quad = null;

    public Light sunLight;

    [SerializeField][HideInInspector]
    public Shader skyShader;
    private Material skyMaterial;

    public float altitude = 1.0f;

    private void Awake()
    {
        if (!quad)
        {
            quad = new Mesh();

            quad.vertices = new Vector3[]{
                new Vector3(-1.0f, -1.0f, 1.0f),
                new Vector3(1.0f, -1.0f, 1.0f),
                new Vector3(-1.0f, 1.0f, 1.0f),
                new Vector3(1.0f, 1.0f, 1.0f)};

            quad.triangles = new int[]{
                0, 2, 1,
                2, 3, 1};

            quad.UploadMeshData(true);
        }
    }

    private void Start()
    {
        skyMaterial = new Material(skyShader);

        atmosphereModel.ComputeLookupTextures();
    }

    DateTime currentTime = DateTime.Now - TimeSpan.FromHours(5); // Time in EST
    double latitude = 45.5017 / (Math.PI * 2.0);
    double longitude = 73.5673 / (Math.PI * 2.0);
    Vector3 sunDirection;

    private void Update()
    {
        // Slows down at sunrise/sunset, and speeds up overnight.
        float sunSpeed = Mathf.Abs(sunDirection.y + 0.1f) * (sunDirection.y < -0.025f ? 10.0f : 1.0f) + 0.2f;
        sunSpeed = 0;
        currentTime = currentTime.AddHours(Time.deltaTime * sunSpeed);
        double sunAltitude;
        double sunAzimuth;
        SunPosition.CalculateSunPosition(currentTime, latitude, longitude, out sunAltitude, out sunAzimuth);
        CartesianCoords.SphericalToCartesian(1, (float)sunAzimuth, (float)sunAltitude, out sunDirection);
        sunLight.transform.rotation = Quaternion.LookRotation(-sunDirection);
        //Vector3 sunDirection = -sunLight.transform.forward;
    }

    private void OnRenderObject()
    {
        const int pass = 0;

        atmosphereModel.SetAtmosphereUniforms(skyMaterial);
        skyMaterial.SetTexture("transmittance_texture", atmosphereModel.TransmittanceLUT);
        skyMaterial.SetTexture("scattering_texture", atmosphereModel.ScatteringLUT);
        skyMaterial.SetTexture("irradiance_texture", atmosphereModel.IrradianceLUT);

        
        Color sunRadiance = new Color(1.0f, 1.0f, 1.0f);
        Vector3 sunSize = new Vector3(1.0f, 1.0f, 1.0f);

        skyMaterial.SetVector("camera", new Vector3(0, (float)(atmosphereModel.PlanetaryRadius / 1000) + altitude, 0));
        skyMaterial.SetVector("sun_direction", sunDirection.normalized);

        if (skyMaterial.SetPass(pass))
        {
            // Create a matrix that will transform verts into camera view directions.
            Matrix4x4 matrix = (GL.GetGPUProjectionMatrix(Camera.current.projectionMatrix, false) * (Camera.current.worldToCameraMatrix * Matrix4x4.TRS(Camera.current.transform.position, Quaternion.identity, Vector3.one))).inverse;
            // Draw the fullscreen quad.
            Graphics.DrawMeshNow(quad, matrix);
        }
        else
            Debug.LogWarning("PhysicalSky: Could not set material pass " + pass + " for rendering");
    }
}
