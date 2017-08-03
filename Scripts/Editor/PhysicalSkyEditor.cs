using UnityEngine;
using UnityEditor;

namespace PhysicalSky
{
    [CustomEditor(typeof(PhysicalSky))]
    public class PhysicalSkyEditor : Editor
    {
        // Toggles are declared static, so they don't keep collapsing between entering/exiting the editor.
        static bool shadersUnfolded = false;

        public override void OnInspectorGUI()
        {
            PhysicalSky physicalSky = serializedObject.targetObject as PhysicalSky;

            if (!physicalSky)
                return;

            physicalSky.SunLightBrightnessMultiplier = EditorGUILayout.FloatField("Sun Light Brightness", physicalSky.SunLightBrightnessMultiplier);
            physicalSky.SunBrightnessMultiplier = EditorGUILayout.FloatField("Sun Brightness", physicalSky.SunBrightnessMultiplier);
            physicalSky.SkyExposure = EditorGUILayout.FloatField("Sky Exposure", physicalSky.SkyExposure);
            physicalSky.Altitude = EditorGUILayout.FloatField("Altitude (km)", physicalSky.Altitude);
            physicalSky.Atmosphere = EditorGUILayout.ObjectField("Atmosphere Model", physicalSky.Atmosphere as AtmosphereModel, typeof(AtmosphereModel), false) as IAtmosphereModel;

            physicalSky.StarBrightnessMultiplier = EditorGUILayout.FloatField("Star Brightness", physicalSky.StarBrightnessMultiplier);
            physicalSky.StarBrightnessPower = EditorGUILayout.FloatField("Star Brightness Power", physicalSky.StarBrightnessPower);
            physicalSky.StarMap = EditorGUILayout.ObjectField("Star Map", physicalSky.StarMap, typeof(StarMap), false) as StarMap;

            shadersUnfolded = EditorGUILayout.Foldout(shadersUnfolded, "Advanced", true);
            if (shadersUnfolded)
            {
                physicalSky.SkyShader = EditorGUILayout.ObjectField("Sky Shader", physicalSky.SkyShader as Shader, typeof(Shader), false) as Shader;
                physicalSky.SunRadianceShader = EditorGUILayout.ObjectField("Sun Radiance Shader", physicalSky.SunRadianceShader as Shader, typeof(Shader), false) as Shader;
                physicalSky.StarMeshShader = EditorGUILayout.ObjectField("Star Mesh Shader", physicalSky.StarMeshShader as Shader, typeof(Shader), false) as Shader;
            }
            EditorUtility.SetDirty(physicalSky);
        }
    }
}