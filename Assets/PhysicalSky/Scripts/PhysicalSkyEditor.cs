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

            physicalSky.Altitude = Mathf.Max(0.001f, EditorGUILayout.FloatField("Altitude (km)", physicalSky.Altitude));
        }
    }
}