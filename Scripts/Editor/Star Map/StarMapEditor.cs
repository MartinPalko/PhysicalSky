using UnityEngine;
using UnityEditor;

namespace PhysicalSky
{
    [CustomEditor(typeof(StarMap))]
    public class StarMapEditor : Editor
    {
        // Toggles are declared static, so they don't keep collapsing between entering/exiting the editor.
        static bool advancedUnfolded = false;

        public override void OnInspectorGUI()
        {
            StarMap starMap = serializedObject.targetObject as StarMap;

            if (!starMap)
                return;

            EditorGUILayout.LabelField("Star map containing " + starMap.Count() + " stars.");

            if (GUILayout.Button("Import from .csv"))
            {
                StarMapImporter mapImporter = EditorWindow.GetWindow<StarMapImporter>(true, "Import star map from CSV");
                mapImporter.starMap = starMap;
                mapImporter.ShowUtility();                
            }

            starMap.FirstIsSun = EditorGUILayout.Toggle("Brightest Star is Sun", starMap.FirstIsSun);
            starMap.BackgroundCube = EditorGUILayout.ObjectField("Background Cubemap", starMap.BackgroundCube, typeof(Cubemap), false) as Cubemap;
            starMap.BackgroundRotation = EditorGUILayout.Vector3Field("Cubemap Rotation", starMap.BackgroundRotation);
            starMap.BackgroundCubeBrightness = EditorGUILayout.FloatField("Cubemap Brightness", starMap.BackgroundCubeBrightness);

            EditorUtility.SetDirty(starMap);
        }
    }
}