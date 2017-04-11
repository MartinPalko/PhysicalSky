using UnityEngine;
using UnityEditor;
using System.IO;
using System.Collections.Generic;

namespace PhysicalSky
{
    public class StarMapImporter : EditorWindow
    {
        public StarMap starMap = null;

        string csvPath = "";
        int numberToImport = 4096;

        class HeaderInfo
        {
            public const string labelX = "x";
            public const string labelY = "y";
            public const string labelZ = "z";
            public const string labelMag = "mag";
            public const string labelCi = "ci";

            public int iX = -1;
            public int iY = -1;
            public int iZ = -1;
            public int iMag = -1;
            public int iCi = -1;

            public HeaderInfo(string[] headerValues)
            {
                for (int i = 0; i < headerValues.Length; i++)
                {
                    string s = headerValues[i];

                    if (s == labelX)
                        iX = i;
                    else if (s == labelY)
                        iY = i;
                    else if (s == labelZ)
                        iZ = i;
                    else if (s == labelMag)
                        iMag = i;
                    else if (s == labelCi)
                        iCi = i;
                }
            }
        }

        void OnGUI()
        {
            if (!starMap)
                return;

            GUILayout.BeginHorizontal();
            csvPath = EditorGUILayout.TextField("Path:", csvPath);
            if (GUILayout.Button("...", GUILayout.Width(30)))
            {
                csvPath = EditorUtility.OpenFilePanelWithFilters("Import star database", Application.dataPath, new string[]{ "Comma seperated values", "csv" });
            }
            GUILayout.EndHorizontal();
            EditorGUILayout.IntField("Number to Import:", numberToImport);

            GUILayout.Space(5);
            if (GUILayout.Button("Import"))
            {
                ImportCSV();
                EditorUtility.SetDirty(starMap);
                Close();
            }
        }

        void ImportCSV()
        {
            StreamReader file = new StreamReader(csvPath);

            string line;
            HeaderInfo headerInfo = null;
            List<StarMap.Star> stars = new List<StarMap.Star>(4096);

            while ((line = file.ReadLine()) != null)
            {
                string[] splitLine = line.Split(',');

                if (headerInfo != null)
                {
                    StarMap.Star newStar = new StarMap.Star();
                    newStar.position = new Vector3(float.Parse(splitLine[headerInfo.iX]), float.Parse(splitLine[headerInfo.iY]), float.Parse(splitLine[headerInfo.iZ]));
                    newStar.colorIndex = float.Parse(splitLine[headerInfo.iCi]);
                    newStar.apparentMagnitude = float.Parse(splitLine[headerInfo.iMag]);
                    stars.Add(newStar);
                }
                else
                    headerInfo = new HeaderInfo(splitLine); // Not set yet, so first line is header info
            }

            // Sort by apparent magnitude
            stars.Sort((a, b) => a.apparentMagnitude.CompareTo(b.apparentMagnitude));
            stars.RemoveRange(numberToImport, stars.Count - numberToImport);

            starMap.Clear();
            starMap.Add(stars);
        }
    }

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
        }
    }
}