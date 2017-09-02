using UnityEngine;
using UnityEditor;
using System.Collections.Generic;
using System.IO;

namespace PhysicalSky
{
    public class StarMapImporter : EditorWindow
    {
        public StarMap starMap = null;

        string csvPath = "";
        int numberToImport = 4096;
        bool zUp = true;
        bool warnOnParseError = false;

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
                csvPath = EditorUtility.OpenFilePanelWithFilters("Import star database", Application.dataPath, new string[] { "Comma seperated values", "csv" });
            }
            GUILayout.EndHorizontal();
            numberToImport = EditorGUILayout.IntField("Number to Import:", numberToImport);
            zUp = EditorGUILayout.Toggle("Z Up:", zUp);

            warnOnParseError = EditorGUILayout.Toggle("Warn on parse error: ", warnOnParseError);

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
                    try
                    {
                        float parsedX = float.Parse(splitLine[headerInfo.iX]);
                        float parsedY = float.Parse(splitLine[headerInfo.iY]);
                        float parsedZ = float.Parse(splitLine[headerInfo.iZ]);

                        newStar.coordinate = new CelestialCoordinates.CartesianCoords(parsedX, zUp ? parsedZ : parsedY, zUp ? parsedY : parsedZ);
                    }
                    catch (System.FormatException e)
                    {
                        Debug.LogWarning("Error parsing star position: " + e.Message);
                        continue;
                    }

                    try
                    {
                        newStar.colorIndex = float.Parse(splitLine[headerInfo.iCi]);
                    }
                    catch (System.FormatException e)
                    {
                        Debug.LogWarning("Error parsing star color index value \'" + splitLine[headerInfo.iCi] + "\': " + e.Message);
                        continue;
                    }

                    try
                    {
                        newStar.apparentMagnitude = float.Parse(splitLine[headerInfo.iMag]);
                    }
                    catch (System.FormatException e)
                    {
                        Debug.LogWarning("Error parsing star color index value \'" + splitLine[headerInfo.iMag] + "\': " + e.Message);
                        continue;
                    }

                    stars.Add(newStar);
                }
                else
                    headerInfo = new HeaderInfo(splitLine); // Not set yet, so first line is header info
            }

            // Sort by apparent magnitude
            stars.Sort((a, b) => a.apparentMagnitude.CompareTo(b.apparentMagnitude));

            int numberToRemove = stars.Count - numberToImport;
            if (numberToRemove > 0)
                stars.RemoveRange(numberToImport, numberToRemove);

            starMap.Clear();
            starMap.Add(stars);
        }
    }
}