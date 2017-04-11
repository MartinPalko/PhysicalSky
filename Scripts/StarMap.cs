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

        [SerializeField]
        private List<Star> m_stars = new List<Star>();

        public void Add(Star star) { m_stars.Add(star); }
        public void Add(IEnumerable<Star> stars) { m_stars.AddRange(stars); }
        public void RemoveAt(int index) { m_stars.RemoveAt(index); }
        public void Clear() { m_stars.Clear(); }
        public int Count() { return m_stars.Count; }
        public Star Get(int index) { return m_stars[index]; }

        public Mesh CreateStarMesh()
        {
            int numVerts = m_stars.Count - (m_stars.Count % 3); // Make sure number of verts is a multiple of 3, since they'll store in triangles.
            Vector3[] verts = new Vector3[numVerts];
            Color[] colors = new Color[numVerts];
            int[] triangles = new int[numVerts];

            for (int i = 0; i < numVerts; i++)
            {
                Star star = m_stars[i];
                verts[i] = star.position;
                // TODO
                //colors[i] = star.color;
                triangles[i] = i; // Just using mesh as a source of verticies to feed the geometry shader, so triangles don't matter.
            }

            Mesh mesh = new Mesh();
            mesh.vertices = verts;
            mesh.colors = colors;
            mesh.triangles = triangles;
            // Make mesh bounds infinitely large, so it is never culled.
            mesh.bounds = new Bounds(Vector3.zero, new Vector3(float.PositiveInfinity, float.PositiveInfinity, float.PositiveInfinity));

            mesh.UploadMeshData(true);

            return mesh;
        }
    }
}