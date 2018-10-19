using UnityEngine;
using UnityEngine.Rendering.PostProcessing;

namespace PhysicalSky
{
    [System.Serializable]
    [PostProcess(typeof(AtmosphericFogRenderer), PostProcessEvent.BeforeTransparent, "Atmospheric Fog", true)]
    public sealed class AtmosphericFog : PostProcessEffectSettings
    {
        public PhysicalSky sky = null;

        public override bool IsEnabledAndSupported(PostProcessRenderContext context)
        {
            sky = FindObjectOfType<PhysicalSky>(); // TEMP HACK

            return enabled.value && sky && sky.enabled && sky.Atmosphere.HasComputedData();
        }
    }

    public sealed class AtmosphericFogRenderer : PostProcessEffectRenderer<AtmosphericFog>
    {
        public override void Render(PostProcessRenderContext context)
        {
            var sheet = context.propertySheets.Get(Shader.Find("Hidden/AtmosphericFogRenderer"));

            settings.sky.SetShaderUniforms(sheet.properties);
            

            var p = GL.GetGPUProjectionMatrix(context.camera.projectionMatrix, true);// Unity flips its 'Y' vector depending on if its in VR, Editor view or game view etc... (facepalm)
            p[2, 3] = p[3, 2] = 0.0f;
            p[3, 3] = 1.0f;
            var clipToWorld = Matrix4x4.Inverse(p * context.camera.worldToCameraMatrix) * Matrix4x4.TRS(new Vector3(0, 0, -p[2, 2]), Quaternion.identity, Vector3.one);
            sheet.properties.SetMatrix("_ClipToWorld", clipToWorld);

            context.command.BlitFullscreenTriangle(context.source, context.destination, sheet, 0);
        }
    }
}