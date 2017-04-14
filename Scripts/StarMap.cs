using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewStarMap", menuName = "PhysicalSky/StarMap")]
    public class StarMap : ScriptableObject
    {
        [System.Serializable]
        public struct Star
        {
            public Vector3 position;
            public float apparentMagnitude;
            public float colorIndex;
        }

        const float SUN_APPARENT_MAGNITUDE = -26.7f;

        [SerializeField]
        private List<Star> m_stars = new List<Star>();

        [SerializeField]
        private bool m_firstIsSun = true;

        public void Add(Star star) { m_stars.Add(star); }
        public void Add(IEnumerable<Star> stars) { m_stars.AddRange(stars); }
        public void RemoveAt(int index) { m_stars.RemoveAt(index); }
        public void Clear() { m_stars.Clear(); }
        public int Count() { return m_stars.Count; }
        public Star Get(int index) { return m_stars[index]; }

        static Color ColorIndexAndMagnitudeToRGB(float colorIndex, float apparentMagnitude)
        {
            float t = 0, r = 0, g = 0, b = 0;

            float bv = Mathf.Clamp(colorIndex, -0.4f, 1.99f);

                if ((bv >= -0.40) && (bv < 0.00)) { t = (bv + 0.40f) / (0.00f + 0.40f); r = 0.61f + (0.11f * t) + (0.1f * t * t); }
            else if ((bv >= 0.00) && (bv < 0.40)) { t = (bv - 0.00f) / (0.40f - 0.00f); r = 0.83f + (0.17f * t); }
            else if ((bv >= 0.40) && (bv < 2.10)) { t = (bv - 0.40f) / (2.10f - 0.40f); r = 1.00f; }

                if ((bv >= -0.40) && (bv < 0.00)) { t = (bv + 0.40f) / (0.00f + 0.40f); g = 0.70f + (0.07f * t) + (0.1f * t * t); }
            else if ((bv >= 0.00) && (bv < 0.40)) { t = (bv - 0.00f) / (0.40f - 0.00f); g = 0.87f + (0.11f * t); }
            else if ((bv >= 0.40) && (bv < 1.60)) { t = (bv - 0.40f) / (1.60f - 0.40f); g = 0.98f - (0.16f * t); }
            else if ((bv >= 1.60) && (bv < 2.00)) { t = (bv - 1.60f) / (2.00f - 1.60f); g = 0.82f - (0.5f * t * t); }

                if ((bv >= -0.40) && (bv < 0.40)) { t = (bv + 0.40f) / (0.40f + 0.40f); b = 1.00f; }
            else if ((bv >= 0.40) && (bv < 1.50)) { t = (bv - 0.40f) / (1.50f - 0.40f); b = 1.00f - (0.47f * t) + (0.1f * t * t); }
            else if ((bv >= 1.50) && (bv < 1.94)) { t = (bv - 1.50f) / (1.94f - 1.50f); b = 0.63f - (0.6f * t * t); }

            float luminosity = Mathf.Pow(10, -apparentMagnitude) / 2.5f;

            return new Color(r, g, b, luminosity);
        }

        public Mesh CreateStarMesh()
        {
            int startIndex = m_firstIsSun ? 1 : 0; // Skip the sun if it's the first in the star list as it's usually way too bright.
            int numVerts = m_stars.Count - startIndex;
            numVerts -= numVerts % 3; // Make sure number of verts is a multiple of 3, since they'll store in triangles.

            Vector3[] verts = new Vector3[numVerts];
            Color[] colors = new Color[numVerts];
            int[] triangles = new int[numVerts];

            for (int i = 0; i < numVerts; i++)
            {
                Star star = m_stars[i + startIndex];
                verts[i] = star.position;
                colors[i] = ColorIndexAndMagnitudeToRGB(star.colorIndex, star.apparentMagnitude);
                triangles[i] = i; // Just using mesh as a source of verticies to feed the geometry shader, so triangles don't matter.
            }

            Mesh mesh = new Mesh();
            mesh.vertices = verts;
            mesh.colors = colors;
            mesh.triangles = triangles;
            // Make mesh bounds extremely large so its not culled.
            mesh.bounds = new Bounds(Vector3.zero, new Vector3(100000.0f, 100000.0f, 100000.0f));

            mesh.UploadMeshData(true);

            return mesh;
        }
    }
}