using UnityEngine;
using UnityEditor;

namespace PhysicalSky
{
    [CustomEditor(typeof(PhysicalSky))]
    public class PhysicalSkyEditor : Editor
    {
        public override void OnInspectorGUI()
        {
            PhysicalSky physicalSky = serializedObject.targetObject as PhysicalSky;

            if (!physicalSky)
                return;

            physicalSky.SunBrightnessMultiplier = EditorGUILayout.FloatField("Sun Light Brightness", physicalSky.SunBrightnessMultiplier);
            physicalSky.SkyExposure = EditorGUILayout.FloatField("Sky Exposure", physicalSky.SkyExposure);
            physicalSky.Altitude = EditorGUILayout.FloatField("Altitude (km)", physicalSky.Altitude);
            physicalSky.Atmosphere = EditorGUILayout.ObjectField("Atmosphere Model", physicalSky.Atmosphere as AtmosphereModel, typeof(AtmosphereModel), false) as IAtmosphereModel;

            physicalSky.StarBrightnessMultiplier = EditorGUILayout.FloatField("Star Brightness", physicalSky.StarBrightnessMultiplier);
            physicalSky.StarBrightnessPower = EditorGUILayout.FloatField("Star Brightness Power", physicalSky.StarBrightnessPower);
            physicalSky.StarMap = EditorGUILayout.ObjectField("Star Map", physicalSky.StarMap, typeof(StarMap), false) as StarMap;

            EditorUtility.SetDirty(physicalSky);

        }
    }
}