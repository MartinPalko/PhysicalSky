using UnityEditor;
using UnityEngine;

namespace PhysicalSky
{
    [CanEditMultipleObjects]
    [CustomEditor(typeof(AtmosphereModel))]
    public class AtmosphereModelEditor : Editor
    {
        public override void OnInspectorGUI()
        {
            DrawDefaultInspector();

            if (GUILayout.Button("Compute Atmosphere"))
            {
                foreach (var t in targets)
                {
                    var atmosphereModel = t as AtmosphereModel;
                    atmosphereModel.Compute(true);
                }
            }

        }
    }
}