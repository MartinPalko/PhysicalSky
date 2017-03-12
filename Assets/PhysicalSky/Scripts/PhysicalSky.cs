using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PhysicalSky : MonoBehaviour
{
    public AtmosphereModel atmosphereModel;
    public GameObject previewMesh = null;
    
    private static Mesh quad = null;

    [SerializeField][HideInInspector]
    public Shader skyShader;
    private Material skyMaterial;

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

    private void OnRenderObject()
    {
        const int pass = 0;

        atmosphereModel.SetAtmosphereUniforms(skyMaterial);
        skyMaterial.SetTexture("transmittance_texture", atmosphereModel.transmittanceLUT);
        skyMaterial.SetTexture("scattering_texture", atmosphereModel.scatteringLUT);
        skyMaterial.SetTexture("irradiance_texture", atmosphereModel.irradianceLUT);

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
