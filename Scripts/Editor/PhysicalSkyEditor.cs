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

            EditorUtility.SetDirty(physicalSky);

        }
    }
}